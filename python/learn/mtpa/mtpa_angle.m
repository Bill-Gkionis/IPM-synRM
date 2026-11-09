function gamma_deg = mtpa_angle(psi_pm, Ld, Lq, I)
%MTPA_ANGLE  The best current angle of the constant-parameter d-q model [deg], from +q toward -d.
%   gamma_deg = mtpa_angle(psi_pm, Ld, Lq, I): psi_pm in Wb, Ld and Lq in H, I the current's peak in A (a vector
%   works too). It solves dT/dgamma = 0 for T = 1.5 p (psi_pm I cos(g) + 0.5 (Lq - Ld) I^2 sin(2g)):
%       psi_pm sin(g) = (Lq - Ld) I cos(2g),   with cos(2g) = 1 - 2 sin(g)^2,
%   a quadratic in sin(g) whose positive root is taken. The same formula as learn/mtpa/report.py in Python.
%   Project note: MTPA - the most torque per amp, section 2.3.
%
%   Example: our band-free rotor with its 10 A constants (FEMM 5.0, fine mesh).
%       mtpa_angle(16.607e-3, 1.1785e-3, 2.1502e-3, [5 10 20 40])   % 14.8, 23.5, 31.7, 37.5 deg

dL = Lq - Ld;
s  = (-psi_pm + sqrt(psi_pm.^2 + 8 * dL.^2 .* I.^2)) ./ (4 * dL .* I);
gamma_deg = asind(s);
end
