% Spherical Bessel function (first kind)
%        h = besseljsph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	h:    Bessel function
%
% See also besselysph, besselhsph, sft, isft.

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function h = besseljsph(m, x)
h = sqrt(pi./(2*x)).*besselj(m+0.5, x);

