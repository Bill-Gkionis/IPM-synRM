% DEMO_PROXIMITY  The core maths of the proximity note in MATLAB: the skin depth, the skin effect of one wire, and the proximity loss, against d / delta.
%
% What it does: for a 0.63 mm copper wire at 100 degC it sweeps the frequency so that d / delta runs from 0.1 to 5 and plots
%   (left)  R_ac / R_dc of the wire alone (rac_rdc_round_wire), and
%   (right) the proximity loss of the wire in a 10 mT (peak) transverse field per unit length, exact (proximity_ratio_exact times proximity_lowfreq) and the
%           low-frequency formula alone.
% It prints the table the Python note prints (python -m optimization.ac_theory), so the two can be compared line by line in the MATLAB validation pass.
%
% NOT RUN YET: written next to the Python (learn/proximity/*.py), to be run in the validation pass. The functions it calls are in this folder:
%   copper_resistivity.m, skin_depth.m, rac_rdc_round_wire.m, proximity_lowfreq.m, proximity_ratio_exact.m.
%
% Expected (copper at 100 degC, from the Python run):
%   d/delta   R_ac/R_dc   proximity ratio
%   0.1       1.0000001   0.9999993
%   0.5       1.0000814   0.9995526
%   1         1.0013007   0.9928916
%   2         1.0204924   0.8976411
%   3         1.0973282   0.6384099
%   5         1.5050394   0.2067688

rho = copper_resistivity(100);
d   = 0.63e-3;
ratio = [0.1 0.2 0.3 0.5 0.7 1 1.5 2 2.5 3 5];                % d / delta
f = rho ./ (pi * 4e-7 * pi .* (d ./ ratio).^2);               % the frequency that gives each d / delta: delta = d / ratio, f = rho / (pi mu0 delta^2)

skin = rac_rdc_round_wire(d, f, rho);
prox = proximity_ratio_exact(d, f, rho);

fprintf('  d/delta        f [Hz]   R_ac/R_dc   proximity ratio\n');
for k = 1:numel(ratio)
    fprintf('  %6.2f  %12.1f   %9.7f   %9.7f\n', ratio(k), f(k), skin(k), prox(k));
end

ff = logspace(2, 6, 400);                                     % the frequencies of the curves
dd = d ./ skin_depth(ff, rho);
figure('Color', [0.988 0.988 0.984]);
subplot(1, 2, 1);
semilogx(dd, rac_rdc_round_wire(d, ff, rho), 'LineWidth', 2);
xlabel('d / \delta'); ylabel('R_{ac} / R_{dc} of one wire alone'); title('Skin effect'); grid on;
subplot(1, 2, 2);
b = 10e-3;                                                    % peak flux density across the wire [T]
loglog(dd, proximity_lowfreq(d, ff, b, rho) .* proximity_ratio_exact(d, ff, rho), 'LineWidth', 2); hold on;
loglog(dd, proximity_lowfreq(d, ff, b, rho), '--', 'LineWidth', 1.5);
xlabel('d / \delta'); ylabel('proximity loss per unit length [W/m] at 10 mT peak'); title('Proximity effect');
legend('exact', 'low-frequency formula', 'Location', 'northwest'); grid on;
