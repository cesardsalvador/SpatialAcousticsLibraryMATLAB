% Derivative of spherical Bessel function
%        dh = besseljsph(m, x)
% Input:
%	m:    order
%	x:    argument
% Output
%	dh:    derivative of spherical Bessel function

% Cesar D. Salvador
% cdsalv@gmail.com

function dh = dbesseljsph(m, x)
dh = (m./x) .* besseljsph(m, x) - besseljsph(m+1, x);
