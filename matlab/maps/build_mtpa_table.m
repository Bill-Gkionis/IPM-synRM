function tbl = build_mtpa_table(M, I_vec, step_deg)
%BUILD_MTPA_TABLE  The MTPA path of the flux maps: for each current magnitude, the current angle that gives the largest torque.
%
%   tbl = build_mtpa_table(M, I_vec, step_deg)
%
%   Written for the maps of the optimized IPM-assisted synRM (fem_maps.mat, possibly rescaled by scale_maps). Validated in MATLAB R2025b on 2026-11-11 (check_fem_maps.m: every check passes).
%
%   WHAT IT DOES, AND WHY
%   A controller is asked for a torque and has to turn it into the current pair (i_d*, i_q*) that makes it with the least current: MTPA, maximum torque per ampere. With the flux
%   maps that is a table, not a formula. For each current magnitude I (A peak) this function walks the circle of that current in the current angle gamma (from the q axis toward -d:
%   i_d = -I sin(gamma), i_q = I cos(gamma), the angle the Simscape page calls the current advance angle), reads the flux linkages from the maps by LINEAR interpolation (interp2, the
%   same as the block's default interpolation and as the Python side's expected values), computes the torque of the maps,
%       T = 1.5 p (psi_d i_q - psi_q i_d)          (amplitude-invariant d-q, p pole pairs),
%   and takes the angle of the largest one. The step in gamma is step_deg (default 0.1 degree); the angle is the sample's own, with no parabola through the neighbours, so that
%   the answer is exactly what the Python side computed for check_fem_maps.m.
%
%   INPUTS   M        the maps struct (fields id_A, iq_A, psi_d_Wb, psi_q_Wb, pole_pairs, rated_current_A; M.scale if it was rescaled)
%            I_vec    the current magnitudes [A peak], default 31 values from 0 to 3 times the rated current of THIS winding (the grid ends there)
%            step_deg the angle step [deg], default 0.1
%   OUTPUT   tbl      struct of row vectors: I_A, gamma_deg, id_A, iq_A, T_Nm (one entry per current)
%
%   USE as a controller's lookup (a torque request T_req [N.m] -> current references):
%       I_ref  = interp1(tbl.T_Nm, tbl.I_A, T_req);          % the current on the path that makes T_req
%       id_ref = interp1(tbl.T_Nm, tbl.id_A, T_req);   iq_ref = interp1(tbl.T_Nm, tbl.iq_A, T_req);
%   (the torque of the path rises with the current, so interp1 may be used; a request beyond the path's last torque returns NaN: the table ends at 3 times the rated current.)
%
%   ERRORS  The circle of a current must stay inside the maps: with the maps' range (i_d from -3 to +1, |i_q| up to 3 times rated) that holds up to 3 times the rated current.

    if nargin < 3 || isempty(step_deg)
        step_deg = 0.1;
    end
    s = 1;
    if isfield(M, 'scale')
        s = M.scale;
    end
    if nargin < 2 || isempty(I_vec)
        I_vec = linspace(0, 3 * M.rated_current_A / s, 31);
    end
    p = M.pole_pairs;
    idA = M.id_A(:)';                    % the d currents: one row of the maps each (increasing)
    iqA = M.iq_A(:)';                    % the q currents: one column each (increasing)
    assert(isequal(size(M.psi_d_Wb), [numel(idA) numel(iqA)]), 'build_mtpa_table:size', 'psi_d_Wb is not one row per i_d and one column per i_q');
    assert(all(diff(idA) > 0) && all(diff(iqA) > 0), 'build_mtpa_table:axes', 'the current vectors must increase');

    gam = 0:step_deg:90;
    n = numel(I_vec);
    tbl.I_A = I_vec(:)';
    tbl.gamma_deg = zeros(1, n);
    tbl.id_A = zeros(1, n);
    tbl.iq_A = zeros(1, n);
    tbl.T_Nm = zeros(1, n);
    for k = 1:n
        I = I_vec(k);
        idq = -I * sind(gam);
        iqq = I * cosd(gam);
        % interp2(X, Y, V, Xq, Yq): V has one row per Y and one column per X, so X = i_q (the columns) and Y = i_d (the rows)
        pd = interp2(iqA, idA, M.psi_d_Wb, iqq, idq, 'linear');
        pq = interp2(iqA, idA, M.psi_q_Wb, iqq, idq, 'linear');
        assert(~any(isnan(pd)) && ~any(isnan(pq)), 'build_mtpa_table:range', 'the %g A circle leaves the flux maps', I);
        T = 1.5 * p * (pd .* iqq - pq .* idq);
        [Tmax, j] = max(T);
        tbl.gamma_deg(k) = gam(j);
        tbl.id_A(k) = idq(j);
        tbl.iq_A(k) = iqq(j);
        tbl.T_Nm(k) = Tmax;
    end
end
