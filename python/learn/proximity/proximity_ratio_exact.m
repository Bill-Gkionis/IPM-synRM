function r = proximity_ratio_exact(d, f_hz, rho)
%PROXIMITY_RATIO_EXACT  The exact proximity loss of a round wire in a uniform transverse field, over the low-frequency formula.
%   r = proximity_ratio_exact(d, f_hz, rho): d the diameter [m], f_hz the frequency, rho the resistivity [Ohm m]. The wire carries no net current.
%   Inside, A = C I1(gamma r) cos(phi); outside, A = (-B r + D / r) cos(phi); matching A and dA/dr at r = a gives C = -2 B / (gamma I0(gamma a)).
%   The power per unit length is (pi sigma w^2 / 2) |C|^2 times the integral from 0 to a of |I1(gamma r)|^2 r dr, and Lommel's integral gives that in
%   closed form, so with q = I1(gamma a) / I0(gamma a) and gamma = (1 + 1j) / delta:
%       p_exact = (2 pi a B^2 / (mu0^2 sigma)) Im[ gamma conj(q) ]
%   Over the low-frequency formula pi sigma w^2 B^2 d^4 / 128 (proximity_lowfreq) that is
%       r = 16 Im[ gamma conj(q) ] / (a^3 k^4),     k^2 = 2 / delta^2,     a = d / 2.
%   1 at low frequency; 0.9929 at d / delta = 1, 0.8976 at 2, 0.6384 at 3 (copper at 100 degC; the Python note prints them from optimization/ac_theory.py).
%
%   Example (not run yet):
%       rho = copper_resistivity(100);  f = 14462.7;  d = 0.63e-3;     % d / delta = 1
%       proximity_ratio_exact(d, f, rho)                               % 0.9928916

a = d ./ 2;
delta = skin_depth(f_hz, rho);
gamma = (1 + 1i) ./ delta;
z = gamma .* a;
q = besseli(1, z, 1) ./ besseli(0, z, 1);
k2 = 2 ./ delta.^2;
r = 16 .* imag(gamma .* conj(q)) ./ (a.^3 .* k2.^2);
end
