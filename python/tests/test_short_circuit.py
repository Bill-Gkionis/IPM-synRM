"""Tests of optimization.short_circuit: the active short circuit's current and torque from the flux maps.

Each test names the property it holds the code to, so a failure says which piece of the physics broke.
"""
import math

import pytest

from optimization import short_circuit as sc
from optimization.maps import FluxMaps


@pytest.fixture(scope="module")
def fm():
    return FluxMaps.load()


def test_braking_torque_is_the_copper_loss_over_the_speed(fm):
    """At zero terminal voltage all the shaft power goes into the copper: T w_m = -1.5 R |i|^2."""
    for rpm in (500, 3000, 16330):
        r = sc.steady_short_circuit(fm, rpm, **sc.HEADLINE)
        assert r["converged"] and r["inside_maps"]
        assert r["torque"] == pytest.approx(r["torque_from_loss"], rel=1e-9)


def test_high_speed_current_tends_to_the_characteristic_current(fm):
    """At high speed R << w L, so the short-circuit current approaches the d current that cancels the magnets' flux."""
    h = sc.HEADLINE
    i_ch = sc.characteristic_current(fm, h["turns"], h["paths"], h["stack_mm"])
    top = sc.steady_short_circuit(fm, 16330, **h)
    assert i_ch < 0
    assert top["i"] == pytest.approx(abs(i_ch), rel=0.01)
    assert abs(top["i_q"]) < 0.02 * top["i"]                                       # almost all of it on -d


def test_characteristic_current_scales_with_the_turns(fm):
    """Ampere-turns set the field: doubling N / a halves the current that cancels the magnets' flux."""
    a = sc.characteristic_current(fm, 9, 4, 60.0)
    b = sc.characteristic_current(fm, 18, 4, 60.0)
    assert b == pytest.approx(a / 2, rel=1e-6)


def test_conduction_speed_formula():
    """sqrt 3 p w_m psi_m = V_dc: 10 mWb, 4 pole pairs, 48 V -> 48 / (sqrt 3 x 0.01 x 4) rad/s."""
    n = sc.conduction_speed(0.01, 4, 48.0)
    assert n == pytest.approx(48.0 / (math.sqrt(3) * 0.01 * 4) * 60 / (2 * math.pi), rel=1e-12)
