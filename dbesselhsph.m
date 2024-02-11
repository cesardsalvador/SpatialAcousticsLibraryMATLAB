% Derivative of spherical Hankel function
%        dh = dbesselhsph(m, k, x)
% Input:
%	m:    order
%	kind: kind
%	x:    argument
% Output
%	dh:   derivative of spherical Hankel function

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 14, 2020

function dh = dbesselhsph(m, kind, x)
dh = (m./x) .* besselhsph(m, kind, x) - besselhsph(m+1, kind, x);