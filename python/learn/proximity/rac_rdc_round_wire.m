function r = rac_rdc_round_wire(d, f_hz, rho)
%RAC_RDC_ROUND_WIRE  R_ac / R_dc of one round wire alone in space (the skin effect), exact.
%   r = rac_rdc_round_wire(d, f_hz, rho): d the diameter [m], f_hz the frequency, rho the resistivity [Ohm m].
%   Inside the wire the current density is J(r) = J0 I0(gamma r) with gamma = (1 + 1j) / delta (I0 the modified Bessel function of order 0),
%   the total current is I = 2 pi J0 a I1(gamma a) / gamma (the integral of J over the disc: the integral of r I0(gamma r) dr is r I1(gamma r) / gamma) and
%   the impedance per unit length is Z' = E(a) / I with E(a) = J0 I0(gamma a) / sigma; divided by the DC resistance per unit length, 1 / (sigma pi a^2),
%       Z' / R_dc' = (z / 2) I0(z) / I1(z),      z = gamma a = (1 + 1j) d / (2 delta),
%   and the real part of that is the AC resistance over the DC one:
%       R_ac / R_dc = Re[ (z / 2) I0(z) / I1(z) ].
%   1 + (d / delta)^4 / 768 at low frequency, 1.0013 at d / delta = 1, 1.0205 at 2, 1.0973 at 3 (copper at 100 degC; the Python note prints them from
%   optimization/ac_theory.py). The scaled Bessel functions (the third argument 1: I_nu(z) exp(-abs(real(z)))) have the same scale factor in both,
%   which cancels in the ratio, so large z does not overflow.
%
%   Example (not run yet):
%       rho = copper_resistivity(100);  f = 14462.7;  d = 0.63e-3;     % d / delta = 1
%       rac_rdc_round_wire(d, f, rho)                                  % 1.0013007

delta = skin_depth(f_hz, rho);
z = (1 + 1i) .* d ./ (2 .* delta);
r = real((z ./ 2) .* besseli(0, z, 1) ./ besseli(1, z, 1));
end
