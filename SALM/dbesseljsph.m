% Derivative of spherical Bessel function
%        dh = besseljsph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	dh:    derivative of spherical Bessel function

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function dh = dbesseljsph(m, x)
dh = (m./x) .* besseljsph(m, x) - besseljsph(m+1, x);
