% Free-field Green function in the spherical harmonic domain
% array
%
%        gnm = gnmff(n, m, rs, x, k, kind, value, norm, source)
%
% Input
%	n:      order
%   m:      degree
%	rs:     source distance
%	x:      observation point in cartesian coordinates [x y z]
%	k:      wave number
%	kind:   kind of spherical Hankel functions
%   value : 'real' for real-valued (default) or 'complex' for complex-valued.
%   norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
%	source: 'pw' (plane wave) or 'sw' (spherical wave)
%
% Output
%	gnm:    spherical-harmonic coefficients of the free-field Green function

% Cesar D. Salvador
% cdsalv@gmail.com

function gnm = gnmff(n, m, rs, x, k, kind, value, norm, source)
[~, ~, r] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
hnrs = besselhsph(n, mod(kind, 2) + 1, k*rs);
jnr = besseljsph(n, k*r);
y = ynm(n, m, x, value, norm);
if strcmp(source, 'pw')
    gnm = 4 * pi * (1j)^n * jnr .* y;
elseif strcmp(source, 'sw')
    gnm = -4 * pi * 1j * k .* jnr .* hnrs * y;
end
