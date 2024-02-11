% Maximum order of a circular harmonic decomposition.
%    N = chorder(r_mic, r_source, k, error)
%
% Input:
% r_mic     : radius of the measurement surface [m].
% r_source	: radius of the sound source [m].
% k         : wave number (2*pi*f/c).
% err       : desired approximation error.
% alpha     : weighting of offset component (low fequency).
% beta      : weighting of curve component (high frequency).
%
% Output:
% N         : maximum order.
%
% See also sst, isst, isht, sht.

function N = chorder(r_mic, r_source, k, err, alpha, beta, p)
a = r_mic/r_source;
Noffset = log(3/2*(1-a)*err) / log(a) + 1;
Ncurve = k*r_mic+0.5*(3*log(1./err) + ...
         0.5*log(k*r_mic)).^(2/3)*(k*r_mic)^(1/3);
N = round(abs( (alpha*Noffset).^p + (beta*Ncurve).^p ).^(1./p));
