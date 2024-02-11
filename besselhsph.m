% Spherical Hankel function
%        h = besselhsph(m, kind, x)
% Input:
%	m:    order
%	kind: kind (1 or 2)
%	x:    argument
% Output
%	h:    Hankel function
%
% See also besseljsph, besselysph, sft, isft.

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 14, 2020

function h = besselhsph(m, kind, x)
if kind == 1
    h = besseljsph(m, x) + 1j*besselysph(m, x);
elseif kind == 2
    h = besseljsph(m, x) - 1j*besselysph(m, x);
end
