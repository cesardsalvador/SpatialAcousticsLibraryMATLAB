% Spherical Bessel function (first kind)
%        h = besseljsph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	h:    Bessel function
%
% See also besselysph, besselhsph, sft, isft.

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 14, 2020

function h = besseljsph(m, x)
h = sqrt(pi./(2*x)).*besselj(m+0.5, x);

