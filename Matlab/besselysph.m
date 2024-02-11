% Spherical Bessel function (second kind)
%        h = besselysph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	h:    Bessel function
%
% See also besseljsph, besselhsph, sft, isft.

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 14, 2020

function h = besselysph(m, x)
h = sqrt(pi./(2*x)).*bessely(m+0.5, x);
