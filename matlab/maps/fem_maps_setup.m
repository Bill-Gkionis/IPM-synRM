% fem_maps_setup.m: load the FEM maps of the optimized IPM-assisted synRM, scale them to a chosen winding, plot them, build the MTPA table and (optionally) fill Simscape's
% FEM-Parameterized PMSM block.
%
% Validated in MATLAB R2025b on 2026-11-11 (sections 1-4; check_fem_maps.m passes all 14 checks). This script is part of the MATLAB side of the flux maps, prepared from the maps' maths and MathWorks' documentation and checked against the Python side
% (optimization/maps_matlab.py writes fem_maps.mat and the numbers a MATLAB run must reproduce). Section 5, the Simulink block, is not run yet.
% Run check_fem_maps.m first in that pass: it compares this chain's MTPA numbers with Python's.
%
% What the file holds (see optimization/maps_matlab.py for every variable): the flux linkages psi_d and psi_q [Wb], the torque [N.m], the 2-D block's Ld, Lq and PM, the 3-D block's
% Fd, Fq and T over a rotor angle vector, the Steinmetz iron-loss tables and the MTPA table, for the 10-turn placeholder winding of the FEMM model: 10 turns a coil side, one path,
% 120 series turns, a 20 mm stack, currents in A PEAK (amplitude-invariant d-q, d = the magnet axis, so L_d < L_q). FEMM 5.0, fine mesh, the maps averaged over a ripple period.
%
% How to use it: set the winding in section 2 (turns per coil side and parallel paths, from the winding study), run section by section (Ctrl+Enter) or all of it.

%% 1. load the maps
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'fem_maps.mat'));                    % the placeholder winding's maps
fprintf('%s\n', S.description);
fprintf('grid: %d values of i_d (%.2f ... %.2f A) x %d of i_q (%.2f ... %.2f A); rated current %.4f A at the placeholder winding, MTPA angle %.2f deg\n', ...
        numel(S.id_A), min(S.id_A), max(S.id_A), numel(S.iq_A), min(S.iq_A), max(S.iq_A), S.rated_current_A, S.rated_gamma_deg);

%% 2. the winding: edit these two numbers
turns_per_coil_side = 20;                                      % N: turns in each coil side (the FEMM model's placeholder is 10)
parallel_paths      = 1;                                       % a: 1, 2 or 4 (the winding study's allowed values)
M = scale_maps(S, turns_per_coil_side, parallel_paths);
fprintf('winding: %g turns a coil side in %g path(s): %g series turns a phase, scale s = %.4g; rated current of this winding %.4f A\n', ...
        M.turns, M.paths, M.series_turns, M.scale, M.rated_current_A / M.scale);

%% 3. look at the maps (x = i_d, y = i_q; the tables are one row per i_d and one column per i_q, so they are transposed for contourf)
figure('Name', 'FEM flux-linkage maps');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; contourf(M.id_A, M.iq_A, 1e3 * M.psi_d_Wb', 20); colorbar; xlabel('i_d [A]'); ylabel('i_q [A]'); title('\psi_d [mWb]'); axis equal tight;
nexttile; contourf(M.id_A, M.iq_A, 1e3 * M.psi_q_Wb', 20); colorbar; xlabel('i_d [A]'); ylabel('i_q [A]'); title('\psi_q [mWb]'); axis equal tight;
nexttile; contourf(M.id_A, M.iq_A, M.T_Nm', 20); colorbar; xlabel('i_d [A]'); ylabel('i_q [A]'); title('mean torque [N.m] (FEMM air-gap integral)'); axis equal tight;
hold on;

%% 4. the MTPA table of this winding
tbl = build_mtpa_table(M);                                     % 31 currents from 0 to 3 times this winding's rated current, a 0.1 degree walk, linear interpolation
plot(tbl.id_A, tbl.iq_A, 'w-', 'LineWidth', 2);                % the MTPA path on the torque map
fprintf('MTPA at 1, 2 and 3 times rated:\n');
for mult = [1 2 3]
    I = mult * M.rated_current_A / M.scale;
    t = build_mtpa_table(M, I);
    fprintf('   %.2f A: gamma %.2f deg, i_d %.2f A, i_q %.2f A, T %.4f N.m\n', I, t.gamma_deg, t.id_A, t.iq_A, t.T_Nm);
end

% a controller's lookup: the torque request to the current references (interp1 on the path; the torque rises with the current)
T_request = 1.5;                                               % N.m
i_ref = interp1(tbl.T_Nm, tbl.I_A, T_request);
id_ref = interp1(tbl.T_Nm, tbl.id_A, T_request);
iq_ref = interp1(tbl.T_Nm, tbl.iq_A, T_request);
fprintf('a request of %.2f N.m: |i| = %.3f A, i_d* = %.3f A, i_q* = %.3f A\n', T_request, i_ref, id_ref, iq_ref);

%% 5. the block (needs a model with a FEM-Parameterized PMSM block; set its path and the winding's resistance first)
% block = 'plant/FEM-Parameterized PMSM';                      % the block's path in the open model
% Rs = 0.05;                                                   % the winding's phase resistance [Ohm], from the copper-loss study: NOT in the maps
% r = fill_fem_pmsm_block(block, M, 'Option', '3d', 'Rs', Rs, 'DryRun', true);   % first a dry run: what it would set and what it cannot find
% r = fill_fem_pmsm_block(block, M, 'Option', '3d', 'Rs', Rs);                   % then for real
% (the plant's own check that the block reproduces the maps: run it at (i_d, i_q) = (tbl.id_A(11), tbl.iq_A(11)) and compare its torque with tbl.T_Nm(11).)
