% b8_mtpa_map.m: the MATLAB half of the bench's case B8 (bench/b8_matlab.py). Python writes flux_map.mat into a folder of its own, runs this script there with
%     matlab -sd <folder> -logfile matlab.log -singleCompThread -noFigureWindows -batch b8_mtpa_map
% and reads answer.mat and answer.json back. Everything this script reads and writes is in that folder (the current folder): it deletes nothing, calls no other program
% and writes nowhere else. Read bench/b8_matlab.py first: it says why the case exists.
%
% What it does
%   1. Loads flux_map.mat: the flux linkages psi_d and psi_q [Wb] that FEMM computed on a grid of d and q currents (id from -10 to 0 A, iq from 0 to 10 A) for the band-free
%      motor, stored one row per id and one column per iq, with the two axes, the pole pairs, the current (10 A) and the angles to search (gamma, degrees).
%   2. Interpolates the map LINEARLY with interp2 ('linear'): the same method as scipy's RegularGridInterpolator('linear'), which Python uses on the same file. The case holds
%      MATLAB's MTPA angle against Python's, so both have to make the same interpolation: it proves the path Python -> .mat -> MATLAB -> answer -> Python, not the method.
%   3. Walks the circle of 10 A in the current angle gamma (id = -I sin(gamma), iq = I cos(gamma)), computes the torque T = 1.5 p (psi_d iq - psi_q id) at each angle, and
%      takes the angle of the largest torque: the MTPA angle (maximum torque per ampere) at 10 A.
%   4. Does the same on a CONTROL map, psi_d = psi_m + Ld id and psi_q = Lq iq, whose MTPA angle has a closed form. If this script read the file's layout the wrong way (rows
%      for columns), the control would show it at once.
%   5. Draws the flux map and the torque curve with MATLAB's own graphics, in MATLAB's light theme (for this figure only), and saves figure.png with exportgraphics (the figure is
%      not visible: no window opens).
%   6. Writes answer.mat (the numbers) and answer.json (the scalars, and facts about this MATLAB: its version, its preferences folder, its threads).
% An error stops the script and MATLAB exits with a nonzero code (matlab -batch), which is how the case sees that it did not finish.

timer = tic;
S = load('flux_map.mat');
needed = {'id_A', 'iq_A', 'psi_d_Wb', 'psi_q_Wb', 'ctrl_psi_d_Wb', 'ctrl_psi_q_Wb', 'pole_pairs', 'current_A', 'gamma_deg'};
for k = 1:numel(needed)
    assert(isfield(S, needed{k}), 'flux_map.mat has no variable %s', needed{k});
end
idA = S.id_A(:)';                     % the d currents, one row of the map each (increasing)
iqA = S.iq_A(:)';                     % the q currents, one column of the map each (increasing)
p = S.pole_pairs;
I = S.current_A;
gam = S.gamma_deg(:)';                % the angles to search, degrees
assert(isequal(size(S.psi_d_Wb), [numel(idA) numel(iqA)]), 'psi_d_Wb is not one row per id and one column per iq');
assert(isequal(size(S.psi_q_Wb), size(S.psi_d_Wb)) && isequal(size(S.ctrl_psi_d_Wb), size(S.psi_d_Wb)) && isequal(size(S.ctrl_psi_q_Wb), size(S.psi_d_Wb)), 'the maps differ in size');
assert(all(diff(idA) > 0) && all(diff(iqA) > 0), 'the axes must increase');

% the circle of I amperes: the current (id, iq) at each angle
idq = -I * sind(gam);
iqq = I * cosd(gam);

[mtpa_deg, mtpa_Nm, torque_Nm] = mtpa_of(idA, iqA, S.psi_d_Wb, S.psi_q_Wb, idq, iqq, gam, p);
[ctrl_mtpa_deg, ctrl_mtpa_Nm, ctrl_torque_Nm] = mtpa_of(idA, iqA, S.ctrl_psi_d_Wb, S.ctrl_psi_q_Wb, idq, iqq, gam, p);

% ---- the picture: MATLAB's own graphics, saved without a window
fig = figure('Visible', 'off', 'Position', [100 100 1500 540]);
% The picture goes to a white page (the project notes' MATLAB chapter), so this one figure gets MATLAB's light theme. What was seen on 2026-11-06: a figure takes the theme MATLAB runs with
% (under a dark graphics theme, set for one session in scratch preferences, the old script gave grey text and a dark panel on white), and with the light theme this one is legible on white.
% MathWorks, https://www.mathworks.com/help/matlab/ref/theme.html (read 2026-11-06): theme(fig,newtheme) "sets the theme of the specified figure", newtheme is "light", "dark", "auto"
% or a GraphicsTheme, "Introduced in R2025a". It changes that figure only: no MATLAB setting is touched.
if ~isMATLABReleaseOlderThan("R2025a")
    theme(fig, "light");
end
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, {sprintf('MATLAB R%s: the flux map of the band-free motor (FEMM 5.0, coarse mesh, %d x %d points), interpolated linearly (interp2)', version('-release'), numel(idA), numel(iqA)), ...
               sprintf('MTPA at %g A: gamma = %.2f deg, T = %.4f N.m (found by walking the %g A circle in %d steps of %.3f deg)', I, mtpa_deg, mtpa_Nm, I, numel(gam), gam(2) - gam(1))}, 'FontSize', 11);
names = {'psi_d [mWb]', 'psi_q [mWb]'};
maps = {S.psi_d_Wb, S.psi_q_Wb};
for k = 1:2
    ax = nexttile(layout);
    contourf(ax, idA, iqA, 1e3 * maps{k}', 14);        % X = id, Y = iq; the map is (id x iq), so Z = its transpose
    hold(ax, 'on');
    plot(ax, idq, iqq, 'w-', 'LineWidth', 2.2);
    plot(ax, -I * sind(mtpa_deg), I * cosd(mtpa_deg), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'MarkerSize', 9);
    hold(ax, 'off');
    colorbar(ax);
    xlabel(ax, 'id [A]'); ylabel(ax, 'iq [A]'); axis(ax, 'equal'); axis(ax, 'tight');
    title(ax, [names{k} ': white arc = the 10 A circle, red = MTPA'], 'FontSize', 10);
end
ax = nexttile(layout);
plot(ax, gam, torque_Nm, 'b-', 'LineWidth', 1.8); hold(ax, 'on');
plot(ax, mtpa_deg, mtpa_Nm, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
hold(ax, 'off'); grid(ax, 'on');
xlim(ax, [0 90]);
xlabel(ax, 'current angle gamma [deg]'); ylabel(ax, 'torque [N.m]');
title(ax, sprintf('T = 1.5 p (psi_d iq - psi_q id) at %g A', I), 'FontSize', 10);
text(ax, 3, 0.18 * max(torque_Nm), {sprintf('MTPA at %.2f deg', mtpa_deg), sprintf('T = %.4f N.m', mtpa_Nm)}, 'HorizontalAlignment', 'left', 'FontSize', 10);
exportgraphics(fig, 'figure.png', 'Resolution', 130);
close(fig);

% ---- the answer: the numbers (.mat) and the scalars with the facts about this MATLAB (.json)
gamma_deg = gam;
save('answer.mat', 'gamma_deg', 'torque_Nm', 'mtpa_deg', 'mtpa_Nm', 'ctrl_torque_Nm', 'ctrl_mtpa_deg', 'ctrl_mtpa_Nm', '-v7');
try
    threads = maxNumCompThreads;
catch
    threads = -1;
end
info = struct('mtpa_deg', mtpa_deg, 'mtpa_Nm', mtpa_Nm, 'ctrl_mtpa_deg', ctrl_mtpa_deg, 'ctrl_mtpa_Nm', ctrl_mtpa_Nm, 'interpolation', 'linear (interp2)', ...
              'matlab_version', version, 'matlab_release', version('-release'), 'matlabroot', matlabroot, 'prefdir', prefdir, 'pwd', pwd, 'max_comp_threads', threads, ...
              'angles', numel(gam), 'map_id_by_iq', size(S.psi_d_Wb), 'seconds', toc(timer));
fid = fopen('answer.json', 'w');
fwrite(fid, jsonencode(info, 'PrettyPrint', true));
fclose(fid);

function [angle, tmax, torque] = mtpa_of(idA, iqA, psi_d, psi_q, idq, iqq, gam, p)
% The torque on the circle and its maximum. interp2(x, y, V, xq, yq): V has one row per y and one column per x, so x = iq (the columns) and y = id (the rows).
    pd = interp2(iqA, idA, psi_d, iqq, idq, 'linear');
    pq = interp2(iqA, idA, psi_q, iqq, idq, 'linear');
    assert(~any(isnan(pd)) && ~any(isnan(pq)), 'the circle leaves the flux map');
    torque = 1.5 * p * (pd .* iqq - pq .* idq);
    [tmax, k] = max(torque);
    angle = gam(k);
end
