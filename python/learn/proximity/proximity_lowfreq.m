function p = proximity_lowfreq(d, f_hz, b_peak, rho)
%PROXIMITY_LOWFREQ  The low-frequency proximity loss of a round wire in a transverse field, per unit length [W/m].
%   p = proximity_lowfreq(d, f_hz, b_peak, rho): d the wire's diameter [m], f_hz the frequency, b_peak the field's PEAK flux density [T] across the
%   wire, rho the resistivity [Ohm m]. The wire carries no net current; the changing field drives an axial electric field E_z = -y dB/dt across it
%   (y measured from the wire's axis, zero net current), so the power per unit length is the integral of sigma E_z^2 over the disc
%       p = sigma (dB/dt)^2 pi d^4 / 64 = pi sigma w^2 B^2 d^4 / 128      for a sinusoid of peak B (the time average of the square is B^2 w^2 / 2).
%   The formula of the winding study (optimization/winding.py), valid while the wire is thin against the skin depth (the exact result over this one
%   is proximity_ratio_exact: 0.993 at d / delta = 1, 0.898 at 2).
%
%   Example (not run yet): a 0.63 mm wire at 20 kHz in 10 mT peak, copper at 100 degC:
%       proximity_lowfreq(0.63e-3, 20e3, 10e-3, copper_resistivity(100))

sigma = 1 ./ rho;
w = 2 * pi .* f_hz;
p = pi .* sigma .* w.^2 .* b_peak.^2 .* d.^4 ./ 128;
end
