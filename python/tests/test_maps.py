"""Tests of the flux maps and the iron-loss map (optimization/iron_loss.py, maps.py, maps_matlab.py). No FEMM and no solver slot: the pure-Python parts, and the saved maps
when they are there (out\\optimization\\maps\\chosen\\flux_map.npz and iron_map.npz, which maps_report.py writes; the tests that need them are skipped without them).

Run from P:   uv run pytest tests/test_maps.py -q      (about 5 s)

In groups:
  1. The loss model: the operators against closed forms, the fit against the datasheet, the thickness scaling.
  2. The winding scaling rule of the maps and of the loss: psi(i) = s psi_0(s i), T(i) = T_0(s i), the same sums at the same ampere-turns.
  3. The maps' own consistency: a two-sided table, the odd and even parity of an old half table, the d-q torque of the maps against the air-gap torque, the MTPA search against brute force.
  4. The MATLAB side's numbers: the expected MTPA of the scaled winding equals the placeholder's at scaled currents, the .mat file reads back.
"""
import math

import numpy as np
import pytest

from optimization import iron_loss, maps, maps_matlab

HAVE_FLUX = maps.FLUX_MAP.exists()
HAVE_IRON = iron_loss.IRON_MAP.exists()


# ---------------------------------------------------------------------------------------------- 1. the loss model
def test_the_waveform_operators_match_their_closed_forms():
    assert iron_loss.self_test(verbose=False)


def test_the_fit_follows_the_datasheet_and_cross_validates():
    f, B, P = iron_loss.datasheet_points()
    assert len(P) == 112
    fit = iron_loss.fit_bertotti(f, B, P)
    assert fit["rms"] < 0.08 and fit["max"] < 0.18                       # the numbers of 2026-11-09 are 6.5 % and 15.4 %
    assert fit["ke"] == pytest.approx(iron_loss.K_E)                     # the classical coefficient is computed, not fitted
    cv = iron_loss.cross_validate()
    assert cv["interior"]["pass"]


def test_the_classical_coefficient_is_faraday_with_the_datasheet_gauge():
    want = math.pi ** 2 * 1.9e6 * (0.356e-3) ** 2 / (6 * 7650.0)
    assert iron_loss.K_E == pytest.approx(want, rel=1e-12)
    assert iron_loss.C_X == pytest.approx(8.763, abs=1e-3)


def test_only_the_classical_sums_change_with_the_lamination():
    assert iron_loss.eddy_scale(0.3) == pytest.approx((0.3 / 0.356) ** 2, rel=1e-12)
    assert iron_loss.eddy_scale(0.356) == pytest.approx(1.0, rel=1e-12)
    with pytest.raises(ValueError):
        iron_loss.eddy_scale(0.0)


@pytest.mark.skipif(not HAVE_IRON, reason="iron_map.npz is made by maps_report.py loss")
def test_the_loss_of_a_thinner_lamination_is_lower_by_the_classical_part_only():
    imap = iron_loss.IronLossMap()
    i_d, i_q, rpm = -4.5, 11.0, 12000.0
    thick, thin = imap(i_d, i_q, rpm, lamination_mm=0.356), imap(i_d, i_q, rpm, lamination_mm=0.3)
    S = imap.sums(i_d, i_q)
    f = 4 * rpm / 60
    classical = sum(S[r, 1] * f ** 2 for r in range(3))
    assert thick["total"] - thin["total"] == pytest.approx((1 - iron_loss.eddy_scale(0.3)) * classical, rel=1e-9)
    assert thick["total"] == pytest.approx(thick["teeth"] + thick["yoke"] + thick["rotor"])


# ---------------------------------------------------------------------------------------------- 2. the winding scaling rule
@pytest.mark.skipif(not HAVE_FLUX, reason="flux_map.npz is made by maps_report.py flux")
def test_scaling_a_winding_rescales_the_flux_and_leaves_the_torque_alone():
    fm = maps.FluxMaps.load()
    i_d, i_q = -3.1, 7.7
    psi0 = fm.psi_dq(2 * i_d, 2 * i_q)                                  # the placeholder at twice the current
    psi = fm.psi_dq(i_d, i_q, turns=20, paths=1)                         # 20 turns: s = 2
    assert psi[0] == pytest.approx(2 * psi0[0], rel=1e-12) and psi[1] == pytest.approx(2 * psi0[1], rel=1e-12)
    assert fm.torque(i_d, i_q, turns=20) == pytest.approx(fm.torque(2 * i_d, 2 * i_q), rel=1e-12)
    assert fm.psi_dq(i_d, i_q, turns=20, paths=2)[0] == pytest.approx(fm.psi_dq(i_d, i_q)[0], rel=1e-12)       # 20 turns in 2 paths is the placeholder
    L0, L = fm.inductances(2 * i_d, 2 * i_q), fm.inductances(i_d, i_q, turns=20)
    assert L["L_qq"] == pytest.approx(4 * L0["L_qq"], rel=1e-9)         # L goes as s^2
    assert fm.torque(i_d, i_q, stack_mm=40.0) == pytest.approx(2 * fm.torque(i_d, i_q), rel=1e-12)    # twice the stack, twice the torque


@pytest.mark.skipif(not HAVE_IRON, reason="iron_map.npz is made by maps_report.py loss")
def test_the_iron_loss_depends_on_the_ampere_turns_only():
    imap = iron_loss.IronLossMap()
    a = imap(-6.0, 9.0, 5000.0, turns=10)
    b = imap(-3.0, 4.5, 5000.0, turns=20)                                # twice the turns at half the current: the same ampere-turns
    c = imap(-6.0, 9.0, 5000.0, turns=20, paths=2)
    assert b["total"] == pytest.approx(a["total"], rel=1e-12) and c["total"] == pytest.approx(a["total"], rel=1e-12)
    assert imap(-6.0, 9.0, 5000.0, stack_mm=40.0)["total"] == pytest.approx(2 * a["total"], rel=1e-12)
    assert imap(-6.0, -9.0, 5000.0)["total"] == pytest.approx(a["total"], rel=1e-12)      # even in i_q
    with pytest.raises(ValueError):
        imap(-60.0, 9.0, 5000.0)                                         # outside the solved region


# ---------------------------------------------------------------------------------------------- 3. the maps' consistency
@pytest.mark.skipif(not HAVE_FLUX, reason="flux_map.npz is made by maps_report.py flux")
def test_the_table_is_two_sided_and_the_torque_of_the_maps_is_the_air_gap_torque():
    fm = maps.FluxMaps.load()
    assert not fm.mirrored and fm.iq_full[0] < 0 < fm.iq_full[-1]
    rng = np.random.default_rng(1)
    i_d, i_q = rng.uniform(-30, 10, 50), rng.uniform(-30, 30, 50)
    T_age, T_dq = fm.torque(i_d, i_q), fm.torque_dq(i_d, i_q)
    assert np.all(np.abs(T_age - T_dq) <= 0.005 * np.abs(T_age) + 0.01)


def test_an_old_half_table_is_mirrored_with_the_right_parity():
    ids, iqs = np.linspace(-10, 0, 5), np.linspace(0, 10, 5)             # a cubic spline needs more than 3 nodes on each axis
    psi_d = np.array([[1.0, 2, 3, 4, 5]] * 5)
    psi_q = np.array([[0.0, 1, 2, 3, 4]] * 5)
    z = {"id_A": ids, "iq_A": iqs, "rated_current_A": 5.0, "rated_gamma_deg": 20.0, "psi_pm0_Wb": 1.0, "pole_pairs": 4, "depth_mm": 20.0,
         "psi_d": psi_d, "psi_q": psi_q, "T_age": psi_q, "T_wst": psi_q, "T_dq": psi_q, "ripple_pp": psi_d}
    fm = maps.FluxMaps(z)
    assert fm.mirrored
    assert fm.psi_dq(-5.0, -10.0)[0] == pytest.approx(fm.psi_dq(-5.0, 10.0)[0])          # psi_d even in i_q
    assert fm.psi_dq(-5.0, -10.0)[1] == pytest.approx(-fm.psi_dq(-5.0, 10.0)[1])         # psi_q odd
    assert fm.torque(-5.0, 0.0) == pytest.approx(0.0, abs=1e-12)                         # and so zero on the d axis


@pytest.mark.skipif(not HAVE_FLUX, reason="flux_map.npz is made by maps_report.py flux")
def test_the_mtpa_search_finds_the_best_angle_of_a_brute_force_scan():
    fm = maps.FluxMaps.load()
    for I in (0.5 * fm.rated, fm.rated, 2.5 * fm.rated):
        g, d, q, T = fm.mtpa(I)
        gam = np.linspace(0, 90, 18001)
        t = fm.torque_dq(-I * np.sin(np.radians(gam)), I * np.cos(np.radians(gam)))
        assert T == pytest.approx(float(t.max()), rel=1e-6)
        assert abs(g - gam[int(np.argmax(t))]) < 0.1


@pytest.mark.skipif(not HAVE_FLUX, reason="flux_map.npz is made by maps_report.py flux")
def test_out_of_range_currents_raise_unless_clipped():
    fm = maps.FluxMaps.load()
    with pytest.raises(ValueError):
        fm.psi_dq(0.0, 4 * fm.rated)
    assert np.isfinite(fm.psi_dq(0.0, 4 * fm.rated, fill="clip")[0])


# ---------------------------------------------------------------------------------------------- 4. the MATLAB side's numbers
@pytest.mark.skipif(not HAVE_FLUX, reason="flux_map.npz is made by maps_report.py flux")
def test_the_expected_mtpa_of_a_scaled_winding_is_the_placeholders_at_scaled_currents():
    M = maps_matlab.build()
    for k in range(3):
        assert M["expected20_gamma_deg"][k] == pytest.approx(M["expected_gamma_deg"][k], abs=1e-9)
        assert M["expected20_T_Nm"][k] == pytest.approx(M["expected_T_Nm"][k], rel=1e-9)
        assert M["expected20_I_A"][k] == pytest.approx(M["expected_I_A"][k] / 2, rel=1e-12)
    ids, iqs = M["id_A"], M["iq_A"]
    a0 = int(np.argmin(np.abs(ids)))
    assert np.allclose(M["Ld_H"] * ids[:, None] + M["PM_Wb"], M["psi_d_Wb"], atol=1e-15)
    cols = np.abs(iqs) > 0
    assert np.allclose(M["Lq_H"][:, cols] * iqs[cols], M["psi_q_Wb"][:, cols], atol=1e-15)
    assert ids[a0] == 0.0 and M["theta_deg"][-1] == pytest.approx(30.0)
    assert M["Fd_Wb"].shape == (len(ids), len(iqs), 5)


@pytest.mark.skipif(not (HAVE_FLUX and HAVE_IRON and maps_matlab.MAT_FILE.exists()), reason="fem_maps.mat is made by maps_matlab.py")
def test_the_mat_file_reads_back_what_was_built():
    import scipy.io
    M = maps_matlab.build()
    back = scipy.io.loadmat(maps_matlab.MAT_FILE)
    for key in ("id_A", "iq_A", "psi_d_Wb", "psi_q_Wb", "T_Nm", "Fd_Wb", "k_Jr_WperHz2", "mtpa_gamma_deg", "expected_T_Nm"):
        want = np.asarray(M[key], float)
        want = want.reshape(1, -1) if want.ndim == 1 else want
        assert np.array_equal(back[key], want), key
