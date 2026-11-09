function delta = skin_depth(f_hz, rho)
%SKIN_DEPTH  The skin depth of a non-magnetic conductor [m]: where an alternating current has fallen to 1/e.
%   delta = skin_depth(f_hz, rho): f_hz the frequency in Hz, rho the resistivity in Ohm m (arrays of one size work).
%   From the diffusion equation of the field in a conductor, laplacian(A) = j w mu0 sigma A, whose solutions decay as exp(-x / delta):
%       delta = sqrt(2 / (w mu0 sigma)) = sqrt(rho / (pi f mu0)),   w = 2 pi f,   sigma = 1 / rho.
%   Same formula as optimization/ac_theory.py (skin_depth) in Python. Copper at 100 degC: rho = copper_resistivity(100) = 2.2662e-8 Ohm m,
%   so delta is 5.36 mm at 200 Hz and 0.536 mm at 20 kHz.
%
%   Example (not run yet: the MATLAB checks come in the validation pass):
%       skin_depth([200 20e3], copper_resistivity(100)) * 1e3      % 5.357  0.536  [mm], the numbers the Python note prints

mu0 = 4e-7 * pi;
delta = sqrt(rho ./ (pi .* f_hz .* mu0));
end
