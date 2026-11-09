function rho = copper_resistivity(t_c)
%COPPER_RESISTIVITY  Annealed copper's resistivity at t_c degrees C [Ohm m].
%   rho = copper_resistivity(t_c): rho_20 (1 + alpha_20 (t_c - 20)) with rho_20 = 1.7241e-8 Ohm m (IEC 60028: 1/58 Ohm mm^2/m) and
%   alpha_20 = 0.00393 per K, the line of the materials guide (02-current.md section 4.3) and of optimization/winding.py and ac_theory.py.
%   2.2662e-8 Ohm m at 100 degC.

rho20   = 1.7241e-8;
alpha20 = 0.00393;
rho = rho20 .* (1 + alpha20 .* (t_c - 20));
end
