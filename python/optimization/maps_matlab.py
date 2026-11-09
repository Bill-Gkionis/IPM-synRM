"""optimization.maps_matlab: the maps in the form Simulink's FEM-Parameterized PMSM block takes: matlab\\maps\\fem_maps.mat, read back and checked in Python.
MATLAB is NOT run here and was not run for it (the MATLAB runs come later, in a validation pass).

What it is, and why
    The plant of the control chain is Simscape Electrical's "FEM-Parameterized PMSM" block. Its tables are exactly what phase 1 computes: flux linkages over the current
    plane, the torque, and the Steinmetz iron-loss coefficients of rotor and stator. This file arranges the maps (flux_map.npz, iron_map.npz) in the block's layouts, adds the
    MTPA table and the numbers a MATLAB run must reproduce, writes them with scipy.io.savemat, and reads the file back with scipy.io.loadmat to prove every array arrived
    unchanged. The MATLAB scripts that use the file are in matlab\\maps\\ (fem_maps_setup.m, scale_maps.m, build_mtpa_table.m, fill_fem_pmsm_block.m, check_fem_maps.m).

The block's parameters (MathWorks, "FEM-Parameterized PMSM", https://www.mathworks.com/help/sps/ref/femparameterizedpmsm.html, read 2026-11-09; R2019b "introduced")
    The page lists the DIALOG names below; it gives NO programmatic (set_param) names, so none are used from memory: fill_fem_pmsm_block.m finds each parameter at run time
    from the dialog prompt in the block's own DialogParameters. Names quoted as the page prints them.
    "Modeling option":  "3-D flux linkage data | No thermal port" (used here for the table of flux linkages and torque; the page also lists 2-D, 3-D and 4-D partial-derivative options,
                        and "| Show thermal port" variants of each)
    "Flux linkage data format": "D and Q axes flux linkages as a function of iD, iQ, theta" (default; others: as a function of I, B, theta; A-phase flux linkage as a function of iD, iQ,
                        theta or of I, B, theta). I = peak current magnitude, B = current advance angle (the angle by which the current leads the q axis, tan^-1(-id/iq): the same angle as
                        this project's gamma).
    "Park's convention for tabulated data": default "Q leads D, rotor angle measured from A-phase to D-axis" (the page's Option 1, the block's internal one); the page's d axis is "the
                        direction of the rotor permanent magnet flux" and the rotor angle is zero when it aligns with the A-phase positive flux: the same as motorgeo.FemmModel.
                        Option 1's transform: i_d = (2/3)(i_a cos(N th) + i_b cos(N th - 2 pi/3) + i_c cos(N th + 2 pi/3)), i_q = -(2/3)(i_a sin(N th) + ...): amplitude-invariant, as FEMM's.
    "Number of pole pairs" (4); "Direct-axis current vector, iD" and "Quadrature-axis current vector, iQ" [A] ("must be two-sided: positive and negative values", so the table holds i_q
                        from -3 to +3 times rated, solved by FEMM on both sides); "Rotor angle vector, theta" [deg] ("four or more points and a range from 0 to 120/N degrees" for the D/Q format);
                        "D-axis flux linkage, Fd(iD,iQ,theta)" and "Q-axis flux linkage, Fq(iD,iQ,theta)" [Wb] (the page says flux linkage is flux times the winding's turns; it does not
                        say whether the magnets' flux is included: the block uses the table as the total flux, with no separate magnet term, and so does this file);
                        "Torque matrix, T(iD,iQ,theta)" [N*m] ("the electromagnetic torque applied to the rotor"); "Interpolation method" (Linear default, or Smooth);
                        "Stator resistance per phase, Rs" [Ohm] (NOT in the maps: the winding's own); "Stator zero-sequence inductance, L0" [H] (enabled for delta or neutral-port windings).
    The 2-D option "sinusoidal back EMF, tabulated with d- and q-axis currents": "Direct axis current vector, id", "Quadrature axis current vector, iq", "Ld matrix, Ld(id,iq)" [H],
                        "Lq matrix, Lq(id,iq)" [H], "Permanent magnet flux linkage, PM(id,iq)" [Wb], "Number of pole pairs", "Stator resistance per phase, Rs", "Stator zero-sequence
                        inductance, L0". Its torque is T = 1.5 N (i_q (i_d L_d + PM) - i_d i_q L_q) (the page's equation as flattened in the text I read: the factor 3/(2N) there is
                        a layout artefact of the page's rendering, 3 N / 2 is the amplitude-invariant d-q torque: CHECK BEFORE CITING against the live page); the page does not say whether
                        L_d, L_q are apparent or incremental: the tables here are the apparent ones so that psi_d = L_d i_d + PM and psi_q = L_q i_q reproduce the flux maps at every node.
    "Iron Losses" -> "Model": "Specify tabulated Steinmetz coefficients": P_rotor(f) = k_hr f + k_Jr f^2 + k_er f^1.5 and P_stator(f) = k_hs f + k_Js f^2 + k_es f^1.5, f the electrical
                        frequency = rotor speed x pole pairs; "Rotor hysteresis loss coefficient, k_hr(id,iq)" [W/Hz], "Rotor eddy current loss coefficient, k_Jr(id,iq)" [W/Hz^2],
                        "Rotor excess current loss coefficient, k_er(id,iq)" [W/Hz^1.5], and the same three for the stator (k_hs, k_Js, k_es) (tables over the d and q current vectors:
                        the same (iD, iQ) grid as the flux tables). Another model of the page, "Tabulate with current and speed": "Rotor iron loss table, Pr(id,iq,w)" [W] with
                        "Rotor speed vector, w" [rpm]: the loss in watts at tabulated speeds, not written here (the three coefficients carry the same information at every speed).
    "Mechanical": "Rotor inertia" [kg*m^2], "Rotor damping" [N*m/(rad/s)]; "Temperature Dependence": "Measurement temperature" [K], "Resistance temperature coefficient" [1/K],
                        "Permanent magnet flux temperature coefficient" (the fractional rate of change of the magnet flux with temperature).
    Not on the page (as read): which index order the 3-D tables must have (the defaults are zeros(3,3,7), implying (current 1, current 2, theta): this file's order), how theta is wrapped or
    extrapolated, and the torque sign convention. These are for the validation pass to confirm.

What is in fem_maps.mat (every array in SI; the winding is the 10-turn placeholder: 10 turns a coil side, one path, 120 series turns; 20 mm stack; scale it with scale_maps.m)
    The tables are on a grid FINE_DIV = 4 times finer than FEMM's (65 i_d values x 97 i_q values, a node every rated / 16 = 0.74 A): the block interpolates LINEARLY by default, and a linear
    interpolation of FEMM's own grid (a node every rated / 4) is up to 0.4 % off in psi between the nodes (maps_report.py, 5.). Each fine table is the bicubic spline of the FEMM grid, which passes
    through every FEMM node, so linear interpolation of the fine table stays within a few hundredths of a per cent of the spline. FEMM's own tables are kept under the names ending _raw.
    id_A (1x65), iq_A (1x97)           the current vectors [A peak] (i_d -35.59 ... +11.86, i_q -35.59 ... +35.59: -3 ... +1 and -3 ... +3 times the rated 11.8627 A)
    psi_d_Wb, psi_q_Wb (65x97)         the flux linkages (row = i_d, column = i_q), the mean over the rotor positions of a ripple period
    T_Nm, T_dq_Nm, ripple_pp_Nm        the mean torque of FEMM's air-gap integral, the d-q torque of the flux maps, and the peak-to-peak ripple (65x97)
    id_raw_A (1x17), iq_raw_A (1x25), psi_d_raw_Wb, psi_q_raw_Wb, T_raw_Nm (17x25)   FEMM's own grid and values, nothing interpolated
    Ld_H, Lq_H, PM_Wb (65x97)          the 2-D block option: psi_d = Ld id + PM, psi_q = Lq iq (PM(id, iq) = psi_d(0, iq), Ld the chord along i_d at that i_q, Lq = psi_q / i_q). Exact at every
                                       node except psi_q on the i_q = 0 column, which the form cannot hold: FEMM's psi_q there is not zero (0.05 ... 0.08 mWb, the mirror asymmetry of the winding)
    theta_deg (1x5), Fd_Wb, Fq_Wb, T3_Nm (65x97x5)   the 3-D block option: the position-averaged maps repeated over theta = 0 ... 120/N = 30 mechanical degrees (a ripple-free table)
    k_hr_WperHz, k_Jr_WperHz2, k_er_WperHz1p5, k_hs_WperHz, k_Js_WperHz2, k_es_WperHz1p5 (65x97)   the Steinmetz tables at the CAD's 0.3 mm lamination (the headline), stator = teeth + yoke
    k_Jr_WperHz2_0p356mm, k_Js_WperHz2_0p356mm                    the classical coefficients at the datasheet's 0.356 mm (the other three are the same at either thickness)
    mtpa_I_A, mtpa_gamma_deg, mtpa_id_A, mtpa_iq_A, mtpa_T_Nm (1x31)   the MTPA path of the maps (current 0 ... 3 x rated), the d-q torque of the interpolated flux maps
    expected_*                         what a MATLAB run must give (see below)
    pole_pairs, turns_ref, paths_ref, series_turns_ref, stack_mm, rated_current_A, rated_gamma_deg, psi_pm0_Wb, description

What a MATLAB run must reproduce (check_fem_maps.m compares; the numbers were made here with the same interpolation as MATLAB's interp2 'linear', the block's default)
    expected_I_A, expected_gamma_deg, expected_T_Nm     the MTPA angle and torque at 1, 2 and 3 times the rated current, found by walking the circle in 0.1 degree steps and taking the
                                                        largest torque 1.5 p (psi_d i_q - psi_q i_d) of the linearly interpolated maps (no refinement), for the placeholder winding;
    expected20_*                                        the same for a 20-turn winding in one path (scale s = 2: the currents are half, the angles and torques the same).

Run (from the Python to FEMM folder, after maps_report's flux and loss sections)
    uv run python -m optimization.maps_matlab           writes matlab\\maps\\fem_maps.mat, reads it back, prints the check and the expected numbers
Numbers to expect: the printout (the arrays' shapes, "read back: every array identical", the expected MTPA values); the file is about 0.4 MB.

Self-check (the answers follow).
  1. The 2-D block option needs L_d(id, iq) with psi_d = L_d i_d + PM. At i_d = 0 that is 0 / 0. What is done, and why does it not change the flux at the node?
  2. Why is the table's i_q axis two-sided although the motor runs with i_q > 0?
  3. The scaled winding has 20 turns in one path. Which tables change value, which only their axes?
Answers.
  1. At i_d = 0 the chord is replaced by the slope d psi_d / d i_d of the interpolated map (its limit). It multiplies i_d = 0, so psi_d = PM there whatever L_d is: the node is
     exact, and L_d is continuous across it.
  2. The block asks for it (the vector "must be two-sided"), and a generating or braking motor has i_q < 0.
  3. Flux linkages and PM double (s = 2) and the inductances quadruple; torque, losses and the angle do not change, only their current axes halve.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.io
from scipy.interpolate import RegularGridInterpolator

import paths
from optimization import iron_loss, maps

MATLAB_DIR = paths.P / "matlab" / "maps"
MAT_FILE = MATLAB_DIR / "fem_maps.mat"
THETA_DEG = np.linspace(0.0, 120.0 / 4, 5)          # mechanical degrees, 0 ... 120 / N (N = 4 pole pairs): five slices, the block wants four or more
LAMINATION_HEADLINE_MM = 0.3
STEP_DEG = 0.1
FINE_DIV = 4                                        # the tables are resampled onto a grid with this many nodes per FEMM interval (rated / 16 apart)


def expected_mtpa(id_A, iq_A, psi_d, psi_q, currents, pp=4, step=STEP_DEG):
    """The MTPA angle and torque by the walk the MATLAB script makes: the circle of each current in `step` degree steps (0 ... 90), the maps interpolated LINEARLY, the largest d-q torque
    taken as it stands (no parabola)."""
    fd = RegularGridInterpolator((id_A, iq_A), psi_d, method="linear")
    fq = RegularGridInterpolator((id_A, iq_A), psi_q, method="linear")
    gam = np.arange(0.0, 90.0 + step / 2, step)
    out = []
    for I in currents:
        d, q = -I * np.sin(np.radians(gam)), I * np.cos(np.radians(gam))
        pts = np.column_stack([d, q])
        T = 1.5 * pp * (fd(pts) * q - fq(pts) * d)
        k = int(np.argmax(T))
        out.append((float(gam[k]), float(T[k])))
    return np.array(out)


def build() -> dict:
    fm = maps.FluxMaps.load()
    imap = iron_loss.IronLossMap()
    ids_raw, iqs_raw = fm.id_A, fm.iq_full                                                  # FEMM's own grid: a node every rated / 4
    # the block interpolates LINEARLY by default (and MATLAB's interp2 'linear' is what check_fem_maps.m uses), and a linear interpolation of the raw grid is 0.4 % off in psi between the
    # nodes (maps_report.py, 5.). So the tables are resampled onto a grid FINE_DIV times finer with the bicubic spline (which passes through every FEMM node): the linear interpolation of the
    # fine tables stays within a few hundredths of a per cent of the spline. The raw FEMM tables go in too (names ending _raw).
    ids = np.linspace(ids_raw[0], ids_raw[-1], (len(ids_raw) - 1) * FINE_DIV + 1)
    iqs = np.linspace(iqs_raw[0], iqs_raw[-1], (len(iqs_raw) - 1) * FINE_DIV + 1)
    ids[np.abs(ids) < 1e-9 * fm.rated] = 0.0                                                # the zero axes stay exactly zero
    iqs[np.abs(iqs) < 1e-9 * fm.rated] = 0.0
    A, B = np.meshgrid(ids, iqs, indexing="ij")
    psi_d, psi_q = fm.psi_dq(A, B)
    T_age, T_dq, ripple = fm.torque(A, B), fm.torque_dq(A, B), fm.ripple(A, B)
    a0 = int(np.argmin(np.abs(ids)))
    # the 2-D option: psi_d = Ld id + PM, psi_q = Lq iq at every node (except the i_q = 0 column for psi_q, see the docstring)
    PM = np.repeat(psi_d[a0:a0 + 1, :], len(ids), axis=0)                                   # PM(id, iq) = psi_d(0, iq)
    L = fm.inductances(A, B)                                                                # the incremental values, for the nodes where a chord is 0 / 0
    with np.errstate(divide="ignore", invalid="ignore"):
        Ld = np.where(np.abs(A) > 1e-9, (psi_d - PM) / A, L["L_dd"])
        Lq = np.where(np.abs(B) > 1e-9, psi_q / B, L["L_qq"])
    # the Steinmetz tables on the same grid
    S = imap.sums(A, B)                                                                     # (3 regions, 3 terms, n_id, n_iq); |i_q| inside
    k_cad = iron_loss.eddy_scale(LAMINATION_HEADLINE_MM)
    stator = S[0] + S[1]
    rotor = S[2]
    # the MTPA table
    I_tab = np.linspace(0.0, 3.0 * fm.rated, 31)
    g_t, d_t, q_t, T_t = fm.mtpa(I_tab[1:])
    exp_I = np.array([1.0, 2.0, 3.0]) * fm.rated
    exp = expected_mtpa(ids, iqs, psi_d, psi_q, exp_I)
    s20 = 20 / 10.0
    exp20 = expected_mtpa(ids / s20, iqs / s20, s20 * psi_d, s20 * psi_q, exp_I / s20)       # the 20-turn winding: scaled axes and flux linkages
    n_th = len(THETA_DEG)
    M = {
        "id_A": ids, "iq_A": iqs, "psi_d_Wb": psi_d, "psi_q_Wb": psi_q, "T_Nm": T_age, "T_dq_Nm": T_dq, "ripple_pp_Nm": ripple,
        "id_raw_A": ids_raw, "iq_raw_A": iqs_raw, "psi_d_raw_Wb": fm.grid["psi_d"], "psi_q_raw_Wb": fm.grid["psi_q"], "T_raw_Nm": fm.grid["T_age"],
        "Ld_H": Ld, "Lq_H": Lq, "PM_Wb": PM,
        "theta_deg": THETA_DEG, "Fd_Wb": np.repeat(psi_d[:, :, None], n_th, axis=2), "Fq_Wb": np.repeat(psi_q[:, :, None], n_th, axis=2),
        "T3_Nm": np.repeat(T_age[:, :, None], n_th, axis=2),
        "k_hr_WperHz": rotor[0], "k_Jr_WperHz2": k_cad * rotor[1], "k_er_WperHz1p5": rotor[2],
        "k_hs_WperHz": stator[0], "k_Js_WperHz2": k_cad * stator[1], "k_es_WperHz1p5": stator[2],
        "k_Jr_WperHz2_0p356mm": rotor[1], "k_Js_WperHz2_0p356mm": stator[1],
        "mtpa_I_A": I_tab, "mtpa_gamma_deg": np.concatenate([[0.0], g_t]), "mtpa_id_A": np.concatenate([[0.0], d_t]), "mtpa_iq_A": np.concatenate([[0.0], q_t]),
        "mtpa_T_Nm": np.concatenate([[0.0], T_t]),
        "expected_I_A": exp_I, "expected_gamma_deg": exp[:, 0], "expected_T_Nm": exp[:, 1],
        "expected20_I_A": exp_I / s20, "expected20_gamma_deg": exp20[:, 0], "expected20_T_Nm": exp20[:, 1],
        "pole_pairs": 4.0, "turns_ref": 10.0, "paths_ref": 1.0, "series_turns_ref": 120.0, "stack_mm": 20.0, "rated_current_A": fm.rated, "rated_gamma_deg": fm.rated_gamma,
        "psi_pm0_Wb": fm.psi_pm0, "lamination_mm": LAMINATION_HEADLINE_MM,
        "description": ("FEMM 5.0 (build 1144), fine mesh, 2D 20 mm, optimized IPM-assisted synRM (36 slots, 8 poles), 10 turns a coil side, one path; currents A peak "
                        "(amplitude-invariant d-q, d = magnet axis, Park option 1); flux linkages Wb; torque N.m; position-averaged over a 60 degree electrical period; "
                        "iron loss: Bertotti fitted to AK Steel M-19 29 gauge, classical term scaled to the 0.3 mm lamination; made by optimization/maps_matlab.py"),
    }
    return M


def main() -> int:
    t0 = time.perf_counter()
    M = build()
    MATLAB_DIR.mkdir(parents=True, exist_ok=True)
    scipy.io.savemat(MAT_FILE, M, do_compression=True)
    print(f"wrote {MAT_FILE} ({MAT_FILE.stat().st_size / 1e3:.0f} kB, {len(M)} variables) in {time.perf_counter() - t0:.2f} s")
    back = scipy.io.loadmat(MAT_FILE)
    bad = []
    for key, value in M.items():
        got = back[key]
        if isinstance(value, str):
            same = str(got[0]) == value
        else:
            want = np.atleast_2d(np.asarray(value, float))
            if np.ndim(value) == 1:
                want = np.asarray(value, float).reshape(1, -1)
            same = got.shape == want.shape and np.array_equal(got, want)
        if not same:
            bad.append(key)
    print("read back with scipy.io.loadmat: " + ("every array identical" if not bad else f"DIFFERENT: {bad}"))
    for key in ("id_A", "iq_A", "psi_d_Wb", "Fd_Wb", "T3_Nm", "k_Jr_WperHz2", "mtpa_I_A"):
        print(f"   {key:22s} shape {back[key].shape}")
    # the 2-D option reproduces the flux maps at every node (and the d-q torque), except psi_q on the i_q = 0 column
    chk_d = np.abs(back["Ld_H"] * back["id_A"].reshape(-1, 1) + back["PM_Wb"] - back["psi_d_Wb"]).max()
    cols = np.abs(back["iq_A"].ravel()) > 0
    chk_q = np.abs(back["Lq_H"][:, cols] * back["iq_A"][:, cols] - back["psi_q_Wb"][:, cols]).max()
    held_back = np.abs(back["psi_q_Wb"][:, ~cols]).max()
    print(f"2-D option: psi_d = Ld id + PM within {chk_d:.1e} Wb and psi_q = Lq iq within {chk_q:.1e} Wb at every node with i_q != 0 (rounding only); on the i_q = 0 column the form "
          f"cannot hold FEMM's psi_q, at most {1e3 * held_back:.3f} mWb ({100 * held_back / np.abs(back['psi_q_Wb']).max():.2f} % of the largest psi_q)")
    # the fine tables against the spline they were made from, and against FEMM's own tables: what linear interpolation costs
    fm = maps.FluxMaps.load()
    lin = back["expected_gamma_deg"].ravel(), back["expected_T_Nm"].ravel()
    cub = [fm.mtpa(I)[0] for I in back["expected_I_A"].ravel()], [fm.mtpa(I)[3] for I in back["expected_I_A"].ravel()]
    raw = expected_mtpa(back["id_raw_A"].ravel(), back["iq_raw_A"].ravel(), back["psi_d_raw_Wb"], back["psi_q_raw_Wb"], back["expected_I_A"].ravel())
    print("what a MATLAB run must reproduce (check_fem_maps.m): MTPA by the 0.1 degree walk on the linearly interpolated (fine) maps; the cubic spline's values and linear on FEMM's own grid for comparison")
    for k, (I, g, T) in enumerate(zip(back["expected_I_A"].ravel(), back["expected_gamma_deg"].ravel(), back["expected_T_Nm"].ravel())):
        print(f"   placeholder winding: I = {I:7.3f} A, gamma {g:6.2f} deg, T {T:.6f} N.m   (cubic spline: {cub[0][k]:.2f} deg, {cub[1][k]:.6f} N.m; linear on the raw grid: {raw[k, 0]:.1f} deg, {raw[k, 1]:.6f} N.m)")
    for I, g, T in zip(back["expected20_I_A"].ravel(), back["expected20_gamma_deg"].ravel(), back["expected20_T_Nm"].ravel()):
        print(f"   20 turns, one path:  I = {I:7.3f} A, gamma {g:6.2f} deg, T {T:.6f} N.m")
    return 1 if (bad or chk_d > 1e-12 or chk_q > 1e-12) else 0


if __name__ == "__main__":
    sys.exit(main())
