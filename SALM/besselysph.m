% Spherical Bessel function (second kind)
%        h = besselysph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	h:    Bessel function
%
% See also besseljsph, besselhsph, sft, isft.

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function h = besselysph(m, x)
h = sqrt(pi./(2*x)).*bessely(m+0.5, x);
