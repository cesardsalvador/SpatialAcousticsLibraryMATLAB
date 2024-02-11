% Maximum order of a spherical harmonic decomposition.
%    N = shorderbound(rm, rs, k, error)
%
% Input:
% rm        : radius of the measurement surface [m].
% rs        : radius of the sound source [m].
% k         : wave number (2*pi*f/c).
% err       : desired approximation error.
%
% Output:
% N         : maximum order.
%
% See also sst, isst, isht, sht.

function N = shorderbound(rm, rs, k, err)
%a = rm/rs;
%Noffset = log(3/2*(1-a)*err) / log(a) + 1;
a = rs/rm;
Noffset = 1./log(a) .* log( a^(3/2) ./ ( (a-1).^(3/2) * err ) ) + 1;
Ncurve = k*rm + 0.5 * ( 3*log(1./err) + 0.5*log(k*rm) ).^(2/3) .* ...
                (k*rm).^(1/3);
p = 4;
N = abs( Noffset.^p + Ncurve.^p ).^(1./p);
