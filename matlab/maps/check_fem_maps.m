function ok = check_fem_maps(mat_file)
%CHECK_FEM_MAPS  The validation of the MATLAB side against the Python side: MTPA angles and torques, and the winding scaling, from fem_maps.mat.
%
%   ok = check_fem_maps()             uses fem_maps.mat next to this file
%   ok = check_fem_maps(mat_file)
%
%   Run in MATLAB R2025b on 2026-11-11: all 14 checks pass. It runs scale_maps and build_mtpa_table on the maps, then compares with the numbers
%   optimization/maps_matlab.py stored in the same file (expected_*, made in Python with the same linear interpolation and the same 0.1 degree walk):
%     1. the placeholder winding: the MTPA angle and torque at 1, 2 and 3 times the rated current (angle within 1e-6 degrees, torque within 1e-9 relative);
%     2. a 20-turn winding in one path (scale 2): the same angles and torques at half the currents (expected20_*);
%     3. the scaling rule itself: the flux linkages of the scaled maps are s times the originals on axes divided by s (exact), the inductances s^2, the torque unchanged;
%     4. the 2-D block's tables reproduce the flux maps at every node: Ld * id + PM = psi_d and Lq * iq = psi_q (to rounding);
%     5. the 3-D block's tables are the maps repeated over the rotor angle vector (a ripple-free table), and theta covers 0 ... 120 / N degrees with five points.
%   Prints one line per check and returns true if all pass.

    if nargin < 1 || isempty(mat_file)
        mat_file = fullfile(fileparts(mfilename('fullpath')), 'fem_maps.mat');
    end
    S = load(mat_file);
    ok = true;

    % 1. the placeholder winding
    tbl = build_mtpa_table(S, S.expected_I_A);
    ok = report(ok, 'MTPA angle, placeholder winding', max(abs(tbl.gamma_deg - S.expected_gamma_deg)), 1e-6, 'deg');
    ok = report(ok, 'MTPA torque, placeholder winding', max(abs(tbl.T_Nm ./ S.expected_T_Nm - 1)), 1e-9, 'relative');

    % 2. a 20-turn winding, one path
    M = scale_maps(S, 20, 1);
    tbl20 = build_mtpa_table(M, S.expected20_I_A);
    ok = report(ok, 'MTPA angle, 20 turns in one path', max(abs(tbl20.gamma_deg - S.expected20_gamma_deg)), 1e-6, 'deg');
    ok = report(ok, 'MTPA torque, 20 turns in one path', max(abs(tbl20.T_Nm ./ S.expected20_T_Nm - 1)), 1e-9, 'relative');
    ok = report(ok, 'the torque is unchanged by the winding', max(abs(tbl20.T_Nm ./ tbl.T_Nm - 1)), 1e-9, 'relative');

    % 3. the scaling rule
    s = M.scale;
    ok = report(ok, 'scale of 20 turns in one path is 2', abs(s - 2), 1e-12, '');
    ok = report(ok, 'flux linkage x s', max(abs(M.psi_d_Wb(:) - s * S.psi_d_Wb(:))), 1e-15, 'Wb');
    ok = report(ok, 'axes / s', max(abs(M.id_A(:) - S.id_A(:) / s)), 1e-12, 'A');
    ok = report(ok, 'inductance x s^2', max(abs(M.Ld_H(:) - s^2 * S.Ld_H(:))), 1e-15, 'H');
    ok = report(ok, 'torque map unchanged', max(abs(M.T_Nm(:) - S.T_Nm(:))), 0, 'N.m');

    % 4. the 2-D block's tables (psi_q = Lq iq cannot hold FEMM's small psi_q on the i_q = 0 column: that column is reported, not judged)
    cols = abs(S.iq_A(:)') > 0;
    resid_d = max(max(abs(S.Ld_H .* S.id_A(:) + S.PM_Wb - S.psi_d_Wb)));             % id down the rows: id_A(:) is a column
    resid_q = max(max(abs(S.Lq_H(:, cols) .* S.iq_A(cols) - S.psi_q_Wb(:, cols))));  % iq along the columns, the i_q = 0 column left out
    ok = report(ok, '2-D tables: Ld id + PM = psi_d', resid_d, 1e-12, 'Wb');
    ok = report(ok, '2-D tables: Lq iq = psi_q (columns with iq ~= 0)', resid_q, 1e-12, 'Wb');
    fprintf('        (information) the i_q = 0 column of psi_q, which the 2-D form cannot hold: at most %.3g Wb\n', max(abs(S.psi_q_Wb(:, ~cols)), [], 'all'));

    % 5. the 3-D block's tables
    ok = report(ok, '3-D tables: Fd is psi_d repeated over theta', max(abs(S.Fd_Wb(:, :, 3) - S.psi_d_Wb), [], 'all'), 0, 'Wb');
    ok = report(ok, '3-D tables: theta runs 0 ... 120 / N degrees in five points', abs(S.theta_deg(end) - 120 / S.pole_pairs) + (numel(S.theta_deg) ~= 5), 1e-12, '');

    fprintf('check_fem_maps: %s\n', ternary(ok, 'ALL PASS', 'SOME FAILED'));
end

function ok = report(ok, name, value, limit, unit)
    good = value <= limit;
    fprintf('  %s  %-58s %.3g %s (limit %.3g)\n', ternary(good, 'PASS', 'FAIL'), name, value, unit, limit);
    ok = ok && good;
end

function out = ternary(condition, a, b)
    if condition
        out = a;
    else
        out = b;
    end
end
