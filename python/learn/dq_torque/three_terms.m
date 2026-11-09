function T = three_terms(p, psi_pm, psi_d, psi_q, i_d, i_q)
%THREE_TERMS  The d-q torque 1.5 p (psi_d i_q - psi_q i_d), split into three field x current products [N.m].
%   T = three_terms(p, psi_pm, psi_d, psi_q, i_d, i_q): p the pole pairs, the flux linkages in Wb, the currents
%   in A (peak, amplitude-invariant d-q). The same split as common.three_terms in Python.
%   Project note: Torque in d-q - field, current and saliency, sections 2.7 and 3.5.
%
%   Example: our band-free rotor at its best angle, 10 A (FEMM 5.0, fine mesh).
%       psi_pm = 16.607e-3; Ld = 1.1785e-3; Lq = 2.1502e-3; id = -3.878; iq = 9.217;
%       T = three_terms(4, psi_pm, psi_pm + Ld*id, Lq*iq, id, iq)

k = 1.5 * p;
T.magnets_x_iq      =  k * psi_pm * i_q;             % the magnets' field, pushed by i_q
T.qfield_x_id       = -k * psi_q * i_d;              % the iron's q field (made by i_q), pushed by i_d
T.dfield_of_id_x_iq =  k * (psi_d - psi_pm) * i_q;   % i_d's own d field (a negative i_d weakens d), pushed by i_q
T.total = T.magnets_x_iq + T.qfield_x_id + T.dfield_of_id_x_iq;
end
