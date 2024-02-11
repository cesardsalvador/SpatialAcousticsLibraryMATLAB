% Derivative of Hankel function
%        dh = dbesselh(m, k, x, s)
% Input:
%	m:    order
%	kind: kind
%	x:    argument
%   s:    scale (0: no scale, 1: by exp(-jx) or exp(jx))
% Output
%	dh:   derivative of Hankel function

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function dh = dbesselh(m, kind, x)
dh = (m./x) .* besselh(m, kind, x) - besselh(m+1, kind, x);