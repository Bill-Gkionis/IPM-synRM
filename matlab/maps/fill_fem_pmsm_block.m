function result = fill_fem_pmsm_block(block, M, varargin)
%FILL_FEM_PMSM_BLOCK  Put the FEM maps into a Simscape Electrical "FEM-Parameterized PMSM" block, finding the block's parameter names at run time.
%
%   result = fill_fem_pmsm_block(block, M, 'Option', '3d', 'Rs', 0.05, ...)
%
%   Written for the maps of the optimized IPM-assisted synRM (fem_maps.mat, rescaled with scale_maps for the chosen winding). NOT RUN YET (it needs a Simulink model with the block; the other scripts were validated in MATLAB R2025b on 2026-11-11): the MATLAB side of the flux maps
%   is prepared and checked by reading. CHECK BEFORE CITING: everything below about the block comes from MathWorks' documentation page of the block
%   (https://www.mathworks.com/help/sps/ref/femparameterizedpmsm.html, read 2026-11-09), which prints the DIALOG names of the parameters and NO programmatic names.
%
%   WHY IT SEARCHES FOR THE NAMES
%   set_param needs a block parameter's programmatic name, which the page does not give and which is not guessed here. A masked block lists its parameters in its
%   'DialogParameters' property: a struct whose fields are the programmatic names, each with a 'Prompt' (the dialog's label). This function matches the documented prompts (below)
%   against those, sets each parameter it finds to the NAME of a base-workspace variable that it creates (the usual way to give a block a table), and reports every prompt it
%   could not find or set, so nothing fails silently. 'DryRun', true lists what it would do and touches neither the block nor the workspace.
%
%   THE PARAMETERS (dialog names as the page prints them; units are the block's defaults: A, Wb, N*m, Ohm, H, deg, kg*m^2)
%   Option '3d'  ("3-D flux linkage data | No thermal port", format "D and Q axes flux linkages as a function of iD, iQ, theta"):
%       "Number of pole pairs" = 4                       "Winding type"                       "Interpolation method" (Linear | Smooth)
%       "Direct-axis current vector, iD" = M.id_A        "Quadrature-axis current vector, iQ" = M.iq_A    (two-sided, as the page requires)
%       "Rotor angle vector, theta" = M.theta_deg [deg: 0 ... 120/N, four or more points]
%       "D-axis flux linkage, Fd(iD,iQ,theta)" = M.Fd_Wb   "Q-axis flux linkage, Fq(iD,iQ,theta)" = M.Fq_Wb   "Torque matrix, T(iD,iQ,theta)" = M.T3_Nm   (index order iD, iQ, theta)
%       "Stator resistance per phase, Rs" (the winding's own: not in the maps; given with 'Rs')
%       "Park's convention for tabulated data": left at the page's default, "Q leads D, rotor angle measured from A-phase to D-axis", which is the convention of the maps
%   Option '2d'  ("2-D partial derivative data | No thermal port", parameterization "sinusoidal back EMF, tabulated with d- and q-axis currents"):
%       "Number of pole pairs", "Direct axis current vector, id" = M.id_A, "Quadrature axis current vector, iq" = M.iq_A,
%       "Ld matrix, Ld(id,iq)" = M.Ld_H, "Lq matrix, Lq(id,iq)" = M.Lq_H, "Permanent magnet flux linkage, PM(id,iq)" = M.PM_Wb, "Stator resistance per phase, Rs"
%   Both options, "Iron Losses" -> Model "Specify tabulated Steinmetz coefficients" with the tables over the same current vectors:
%       "Rotor hysteresis loss coefficient, k_hr(id,iq)" [W/Hz] = M.k_hr_WperHz,  "Rotor eddy current loss coefficient, k_Jr(id,iq)" [W/Hz^2] = M.k_Jr_WperHz2,
%       "Rotor excess current loss coefficient, k_er(id,iq)" [W/Hz^1.5] = M.k_er_WperHz1p5, and the stator's k_hs, k_Js, k_es (M.k_hs_WperHz, M.k_Js_WperHz2, M.k_es_WperHz1p5)
%       (the maps give k_J at the CAD's 0.3 mm lamination; the datasheet's 0.356 mm values are M.k_Jr_WperHz2_0p356mm and M.k_Js_WperHz2_0p356mm)
%   Optional, through name-value pairs: 'Inertia' ("Rotor inertia" [kg*m^2]), 'Damping' ("Rotor damping" [N*m/(rad/s)]).
%   The enumerated parameters ("Modeling option", "Flux linkage data format", "Parameterization", the iron-loss "Model") are set from the strings the page prints; if the block's
%   internal enumeration strings differ from the labels, set_param reports it and the entry appears in result.failed: set those four by hand once and rerun.
%
%   INPUTS   block   the path of the block in an open model, e.g. 'plant/FEM-Parameterized PMSM'
%            M       the maps struct (S = load('fem_maps.mat'); M = scale_maps(S, N, a))
%   NAME-VALUE   'Option' '3d' (default) | '2d';  'Rs' [Ohm];  'Interpolation' 'Linear' (default) | 'Smooth';  'WindingType' 'Wye-wound' (default);
%                'Inertia', 'Damping';  'IronLoss' true (default) | false;  'DryRun' false (default) | true
%   OUTPUT   result  struct: set (cell array of {prompt, parameter name, variable}), notfound (prompts), failed ({prompt, parameter name, message}), prompts (all the block's prompts)

    p = inputParser;
    addParameter(p, 'Option', '3d', @(x) any(strcmp(x, {'3d', '2d'})));
    addParameter(p, 'Rs', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'Interpolation', 'Linear', @ischar);
    addParameter(p, 'WindingType', 'Wye-wound', @ischar);
    addParameter(p, 'Inertia', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'Damping', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'IronLoss', true, @islogical);
    addParameter(p, 'DryRun', false, @islogical);
    parse(p, varargin{:});
    o = p.Results;

    dp = get_param(block, 'DialogParameters');
    assert(isstruct(dp) && ~isempty(fieldnames(dp)), 'fill_fem_pmsm_block:dialog', 'the block %s has no DialogParameters: is it the FEM-Parameterized PMSM block?', block);
    names = fieldnames(dp);
    prompts = cell(size(names));
    for k = 1:numel(names)
        entry = dp.(names{k});
        if isstruct(entry) && isfield(entry, 'Prompt')
            prompts{k} = entry.Prompt;
        else
            prompts{k} = names{k};
        end
    end

    result.set = {};
    result.notfound = {};
    result.failed = {};
    result.prompts = prompts;

    pairs = 4;                                                  % number of pole pairs of the maps (M.pole_pairs)
    if isfield(M, 'pole_pairs')
        pairs = M.pole_pairs;
    end

    % each row: documented prompt, value, base-workspace variable name ('' = set the value as the text it is)
    rows = {};
    if strcmp(o.Option, '3d')
        rows(end+1, :) = {'Modeling option', '3-D flux linkage data | No thermal port', ''};
        rows(end+1, :) = {'Flux linkage data format', 'D and Q axes flux linkages as a function of iD, iQ, theta', ''};
        rows(end+1, :) = {'Winding type', o.WindingType, ''};
        rows(end+1, :) = {'Number of pole pairs', pairs, 'fem_pole_pairs'};
        rows(end+1, :) = {'Direct-axis current vector, iD', M.id_A, 'fem_iD'};
        rows(end+1, :) = {'Quadrature-axis current vector, iQ', M.iq_A, 'fem_iQ'};
        rows(end+1, :) = {'Rotor angle vector, theta', M.theta_deg, 'fem_theta'};
        rows(end+1, :) = {'D-axis flux linkage, Fd(iD,iQ,theta)', M.Fd_Wb, 'fem_Fd'};
        rows(end+1, :) = {'Q-axis flux linkage, Fq(iD,iQ,theta)', M.Fq_Wb, 'fem_Fq'};
        rows(end+1, :) = {'Torque matrix, T(iD,iQ,theta)', M.T3_Nm, 'fem_T'};
        rows(end+1, :) = {'Interpolation method', o.Interpolation, ''};
    else
        rows(end+1, :) = {'Modeling option', '2-D partial derivative data | No thermal port', ''};
        rows(end+1, :) = {'Parameterization', 'sinusoidal back EMF, tabulated with d- and q-axis currents', ''};
        rows(end+1, :) = {'Winding type', o.WindingType, ''};
        rows(end+1, :) = {'Number of pole pairs', pairs, 'fem_pole_pairs'};
        rows(end+1, :) = {'Direct axis current vector, id', M.id_A, 'fem_id'};
        rows(end+1, :) = {'Quadrature axis current vector, iq', M.iq_A, 'fem_iq'};
        rows(end+1, :) = {'Ld matrix, Ld(id,iq)', M.Ld_H, 'fem_Ld'};
        rows(end+1, :) = {'Lq matrix, Lq(id,iq)', M.Lq_H, 'fem_Lq'};
        rows(end+1, :) = {'Permanent magnet flux linkage, PM(id,iq)', M.PM_Wb, 'fem_PM'};
        rows(end+1, :) = {'Interpolation method', o.Interpolation, ''};
    end
    if ~isempty(o.Rs)
        rows(end+1, :) = {'Stator resistance per phase, Rs', o.Rs, 'fem_Rs'};
    end
    if o.IronLoss
        rows(end+1, :) = {'Model', 'Specify tabulated Steinmetz coefficients', ''};
        rows(end+1, :) = {'Rotor hysteresis loss coefficient, k_hr(id,iq)', M.k_hr_WperHz, 'fem_khr'};
        rows(end+1, :) = {'Rotor eddy current loss coefficient, k_Jr(id,iq)', M.k_Jr_WperHz2, 'fem_kJr'};
        rows(end+1, :) = {'Rotor excess current loss coefficient, k_er(id,iq)', M.k_er_WperHz1p5, 'fem_ker'};
        rows(end+1, :) = {'Stator hysteresis loss coefficient, k_hs(id,iq)', M.k_hs_WperHz, 'fem_khs'};
        rows(end+1, :) = {'Stator eddy current loss coefficient, k_Js(id,iq)', M.k_Js_WperHz2, 'fem_kJs'};
        rows(end+1, :) = {'Stator excess current loss coefficient, k_es(id,iq)', M.k_es_WperHz1p5, 'fem_kes'};
    end
    if ~isempty(o.Inertia)
        rows(end+1, :) = {'Rotor inertia', o.Inertia, 'fem_J'};
    end
    if ~isempty(o.Damping)
        rows(end+1, :) = {'Rotor damping', o.Damping, 'fem_B'};
    end

    normalise = @(t) lower(regexprep(t, '[^a-zA-Z0-9]', ''));
    keys = cellfun(normalise, prompts, 'UniformOutput', false);
    for r = 1:size(rows, 1)
        wanted = rows{r, 1};
        value = rows{r, 2};
        var = rows{r, 3};
        hit = find(strcmp(keys, normalise(wanted)));
        if isempty(hit)
            result.notfound{end+1} = wanted;
            continue
        end
        for h = hit(:)'
            pname = names{h};
            if o.DryRun
                result.set(end+1, :) = {wanted, pname, var};
                continue
            end
            try
                if isempty(var)
                    set_param(block, pname, value);
                else
                    assignin('base', var, value);
                    set_param(block, pname, var);
                end
                result.set(end+1, :) = {wanted, pname, var};
            catch err
                result.failed(end+1, :) = {wanted, pname, err.message};
            end
        end
    end

    fprintf('fill_fem_pmsm_block (%s, option %s%s): %d parameter(s) set, %d prompt(s) not found, %d failed\n', block, o.Option, ...
            ternary(o.DryRun, ', DRY RUN', ''), size(result.set, 1), numel(result.notfound), size(result.failed, 1));
    for k = 1:numel(result.notfound)
        fprintf('   not found in the block: "%s"\n', result.notfound{k});
    end
    for k = 1:size(result.failed, 1)
        fprintf('   could not set "%s" (%s): %s\n', result.failed{k, 1}, result.failed{k, 2}, result.failed{k, 3});
    end
end

function out = ternary(condition, a, b)
    if condition
        out = a;
    else
        out = b;
    end
end
