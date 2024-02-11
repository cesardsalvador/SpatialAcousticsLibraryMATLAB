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

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% Reference and citation
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.
%     DOI: 10.1109/TASLP.2017.2778562
% [2] C. D. Salvador et al., “Design theory for binaural synthesis:
%     Combining microphone array recordings and head-related transfer function datasets,”
%     Acoust. Sci. Technol., vol. 38, no. 2, pp. 51–62, Mar. 2017.
%     DOI: 10.1250/ast.38.51


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
