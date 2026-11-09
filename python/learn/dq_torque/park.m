function [d, q] = park(xa, xb, xc, theta_e_deg)
%PARK  Phase quantities to the rotor's d-q frame, amplitude-invariant.
%   [d, q] = park(xa, xb, xc, theta_e_deg) takes the three phase values (currents in A, or flux linkages in Wb)
%   and the rotor's electrical angle in degrees, and returns their d and q parts. It is the same transform as
%   motorgeo's FemmModel.solve_abc in Python: a phase peak of 10 A gives sqrt(d^2 + q^2) = 10 A.
%   Project note: Torque in d-q - field, current and saliency, section 3.1.
%
%   Example: the best current at 10 A, through the inverse transform (what FemmModel.solve sends to FEMM) and back.
%       th   = deg2rad(37 + [0 -120 120]);              % any rotor angle: 37 deg electrical here
%       iabc = -3.878 * cos(th) - 9.217 * sin(th);      % i_d cos - i_q sin, phase by phase
%       [id, iq] = park(iabc(1), iabc(2), iabc(3), 37)  % gives -3.878 and 9.217 back

th = deg2rad(theta_e_deg + [0 -120 120]);   % each phase's axis, seen from the rotor
x  = [xa xb xc];
d  =  2/3 * sum(x .* cos(th));              % the part along the magnets (d)
q  = -2/3 * sum(x .* sin(th));              % the part 90 electrical degrees ahead of it (q)
end
