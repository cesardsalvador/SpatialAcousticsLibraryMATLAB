% Sphere transfer function
%   H = sphere_transfer_function(a, b, theta, k, thres)
%
% INPUT
% a:        radius of the measurement point (mic)
% b:        radius of the radiating point (source)
% theta:    incident angle between mic and source
% k:        wave number (k = 2*pi*frequency/velocity_of_sound_in_air)
% thres:    threshold error (typical 0.0001)
%
% OUPUT
% H:        sphere transfer function
%
% See also pressure_rigid_sphere.

% Reference
% [1] R. O. Duda y W. L. Martens, "Range dependence of the
% response of a spherical head model", J. Acoust. Soc. Am.,
% vol. 104, n. 5, pp. 3048-3058, nov. 1998.
% [2]J. J. Bowman, T. B. A. Senior, and P. L. E. Uslenghi,
% Electromagnetic and acoustic scattering by simple shapes,
% New York, NY, USA: Hemisphere, 1987.

% Modified by Cesar D. Salvador
% cesardsalvador@gmail.com

function [H, n] = sphere_transfer_function(a, b, theta, k, thres)

x = cosd(theta);
mu = k * a;
rho = b / a;

zr = 1 ./ (1j*mu * rho);
za = 1 ./ (1j * mu);

Qr2 = zr;
Qr1 = zr .* (1 - zr);
Qa2 = za;
Qa1 = za .* (1 - za);

P2 = 1;
P1 = x;

total = 0;
term = zr / (za * (za - 1));
total = total + term;
term = (3*x * zr * (zr - 1)) / (za *(2 * za^2 - 2 * za +1 ));
total = total + term;

oldratio = 1;
newratio = abs(max(abs(term(:)))) / abs(max(abs(total(:))));

m = 2;

while (oldratio > thres) || (newratio > thres)
    Qr = -(2*m - 1)* zr * Qr1 + Qr2;
    Qa = -(2*m - 1)* za * Qa1 + Qa2;
    P = ( (2*m - 1)* x .* P1 - (m - 1) * P2 )/m;
    term = ( (2*m + 1)* P * Qr ) / ( (m + 1)* za * Qa - Qa1);
    total = total + term;
    m = m + 1;
    Qr2 = Qr1;
    Qr1 = Qr;
    Qa2 = Qa1;
    Qa1 = Qa;
    P2 = P1;
    P1 = P;
    oldratio = newratio;
    newratio = abs(max(abs(term(:)))) / abs(max(abs(total(:))));
end

H = (rho * exp(-1j * mu) * total) ./ (1j * mu);
H = conj(H);
n = m - 1;
