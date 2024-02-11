% Rigid sphere with incident wave (spherical harmonic domain model)
%
%        pnm = rigidspherenm(n, m, a, xs, source, k, kind, value, norm)
%
% Input
%	n:      order
%   m:      degree
%   a:      radius of rigid sphere
%	xs:     source position (only angle is used in case of plane-wave)
%	source: 'pw' (plane wave) or 'sw' (spherical wave)
%	k:      wave number
%	kind:   kind of spherical Hankel functions
%   value : 'real' for real-valued (default) or 'complex' for complex-valued.
%   norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
%
% Output
%	pwnm:   spherical-harmonic coefficients of the free-field Green function

% Cesar D. Salvador
% cdsalv@gmail.com

function pnm = rigidspherenm(n, m, a, xs, source, k, kind, value, norm)
[~, ~, rs] = cart2sph(xs(1, 1), xs(1, 2), xs(1, 3));
dhna = dbesselhsph(n, kind, k*a);
hnrs = besselhsph(n, kind, k*rs);
cys = conj(ynm(n, m, xs, value, norm));
if strcmp(source, 'pw')
    pnm = 4*pi * (1j)^(n+1) * cys ./ ( (k*a).^2 .* dhna );
elseif strcmp(source, 'sw')
    pnm = -4*pi * hnrs * cys ./ ( k * a^2 .* dhna );
end
