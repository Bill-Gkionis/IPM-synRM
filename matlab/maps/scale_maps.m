function M = scale_maps(M0, turns_per_coil_side, parallel_paths)
%SCALE_MAPS  Scale the FEM maps of the placeholder winding to another winding.
%
%   M = scale_maps(M0, turns_per_coil_side, parallel_paths)
%
%   Written for the maps of the optimized IPM-assisted synRM (FEMM 5.0, fine mesh): fem_maps.mat, loaded into the struct M0 with S = load('fem_maps.mat').
%   NOT RUN: the MATLAB side of the flux maps is prepared and checked by reading, not by running (the MATLAB runs come later, in a validation pass).
%
%   WHAT IT DOES, AND WHY
%   FEMM multiplies the current of a winding block by its turns, so the magnetic field depends on the AMPERE-TURNS only. The maps were solved for the placeholder winding: 10
%   turns a coil side, one path, 12 coils a phase in series, so 120 series turns (M0.turns_ref, M0.paths_ref, M0.series_turns_ref). Another winding with N turns a coil side in
%   `a` parallel paths has the series turns 12 N / a, and carries the phase current i. With the scale
%       s = (N / a) / (turns_ref / paths_ref)
%   its current i produces the field of the placeholder's current s*i, and its flux linkage is s times the placeholder's (more series turns link the same flux):
%       i_d' = i_d0 / s,   psi' = s psi0,   L' = s^2 L0,   T' = T0,   losses (watts) unchanged
%   so the SAME tables are valid on rescaled axes: the current vectors are divided by s, the flux linkages (psi_d, psi_q, PM, Fd, Fq) multiplied by s, the inductances (Ld, Lq) by
%   s^2; the torque, the ripple, the iron-loss coefficients and the angles are not touched, and the MTPA table keeps its angles and torques with its currents divided by s. The winding
%   study checked the rule bit for bit in FEMM (5, 10 and 20 turns at 20, 10 and 5 A gave identical torque and field).
%   What the maps do NOT scale: the stator resistance Rs (the winding's own) and the end-winding leakage (a 2D model has none): add them in the plant.
%
%   INPUTS   M0 the struct of fem_maps.mat; turns_per_coil_side N (> 0); parallel_paths a (> 0, default 1)
%   OUTPUT   M  a copy of M0 with the fields above scaled, plus M.scale (= s), M.turns, M.paths and M.series_turns (= series_turns_ref * s)
%
%   EXAMPLE  (placeholder values of s: 20 turns in one path -> s = 2; 20 turns in two paths -> s = 1)
%       S = load('fem_maps.mat');  M = scale_maps(S, 20, 1);   % currents halve, flux linkages double
%       max(abs(M.psi_d_Wb(:) - 2 * S.psi_d_Wb(:)))            % 0

    if nargin < 3 || isempty(parallel_paths)
        parallel_paths = 1;
    end
    assert(isscalar(turns_per_coil_side) && turns_per_coil_side > 0, 'scale_maps:turns', 'turns_per_coil_side must be a positive scalar');
    assert(isscalar(parallel_paths) && parallel_paths > 0, 'scale_maps:paths', 'parallel_paths must be a positive scalar');
    required = {'turns_ref', 'paths_ref', 'series_turns_ref', 'id_A', 'iq_A', 'psi_d_Wb', 'psi_q_Wb'};
    for k = 1:numel(required)
        assert(isfield(M0, required{k}), 'scale_maps:field', 'the maps have no field %s', required{k});
    end

    s = (turns_per_coil_side / parallel_paths) / (M0.turns_ref / M0.paths_ref);
    M = M0;

    % the current axes (and the MTPA table's currents): the winding's current is the placeholder's divided by s
    for name = {'id_A', 'iq_A', 'id_raw_A', 'iq_raw_A', 'mtpa_I_A', 'mtpa_id_A', 'mtpa_iq_A', 'expected_I_A'}
        if isfield(M0, name{1})
            M.(name{1}) = M0.(name{1}) / s;
        end
    end
    % flux linkages go with the series turns
    for name = {'psi_d_Wb', 'psi_q_Wb', 'PM_Wb', 'Fd_Wb', 'Fq_Wb', 'psi_d_raw_Wb', 'psi_q_raw_Wb'}
        if isfield(M0, name{1})
            M.(name{1}) = s * M0.(name{1});
        end
    end
    % inductances go with the turns squared
    for name = {'Ld_H', 'Lq_H'}
        if isfield(M0, name{1})
            M.(name{1}) = s^2 * M0.(name{1});
        end
    end

    M.scale = s;
    M.turns = turns_per_coil_side;
    M.paths = parallel_paths;
    M.series_turns = M0.series_turns_ref * s;
end
