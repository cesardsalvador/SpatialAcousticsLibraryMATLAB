% Derivative of spherical Hankel function
%        dh = dbesselhsph(m, k, x)
% Input:
%	m:    order
%	kind: kind
%	x:    argument
% Output
%	dh:   derivative of spherical Hankel function

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function dh = dbesselhsph(m, kind, x)
dh = (m./x) .* besselhsph(m, kind, x) - besselhsph(m+1, kind, x);