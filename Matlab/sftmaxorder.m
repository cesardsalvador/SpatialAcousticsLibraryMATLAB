% Maximum order of a spherical Fourier transform
%    N = sftmaxorder(a, b, k, error)
%
% Input:
% a        : radius of the measurement surface (m)
% b        : radius of the sound source (m)
% k        : wave number (2*pi*f/c)
% err      : specified order truncation error (e.g. 1e-3)
%
% Output:
% N        : maximum order.
%
% See also sst, isst, isht, sht.

% Cesar D. Salvador
% cdsalv@gmail.com

function N = sftmaxorder(a, b, k, err)
ka = k*a;
r = b/a;
p = 4;
Ncurve = ka + 0.5 * ( 3*log(1./err) + 0.5*log(ka) ).^(2/3) .* (ka).^(1/3);
Noffset = 1./log(r) .* log( r^(3/2) ./ ( (r-1).^(3/2) * err ) ) + 1;
N = floor( abs( Ncurve.^p + Noffset.^p ).^(1./p) );
