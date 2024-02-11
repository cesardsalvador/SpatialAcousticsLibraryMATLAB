% Transfer function for a plane wave impinging on the surface of a rigid sphere.
%   TF = transferFunctionRigidSpherePlaneWave(a, theta, k, thres)
%
% INPUT
%   a       : radius of the point of measurement on the rigid spherical surface
%   theta   : incident angle between the point of measurement and the point source
%   k       : wave number (k = 2*pi*frequency/velocity_of_sound_in_air)
%   thres   : threshold error (typical 0.0001)
%
% OUPUT
%   TF      : Transfer function on the point of measurement due to an incident plane wave
%
% See also: transferFunctionRigidSphere, pressureRigidSphere, pressureRigidSpherePlaneWave

% Reference
% [1] J. J. Bowman, T. B. A. Senior, and P. L. E. Uslenghi,
% Electromagnetic and acoustic scattering by simple shapes,
% chapter 10.3.2, pp. 374-376, New York, NY, USA: Hemisphere, 1987.

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 12, 2020

function TF = transferFunctionRigidSpherePlaneWave(a, theta, k, thres)
kind = 1;                   % kind of spherical Hankel functions (1 or 2)
cosTheta = cosd(theta);     % cosine of incidence angle

n = 0;
total = 0;
term = (2*n+1) * pnm(n, 0, cosTheta) ./ dbesselhsph(n, kind, k*a);
total = total + term;

n = 1;
term = -1j * (2*n+1) * pnm(n, 0, cosTheta) ./ dbesselhsph(n, kind, k*a);
total = total + term;

oldratio = 1;
newratio = abs(max(abs(term(:))))/abs(max(abs(total(:))));

n = 2;
while (oldratio > thres) || (newratio > thres)
    term = (-1j)^n * (2*n+1) * pnm(n, 0, cosTheta) ./ dbesselhsph(n, kind, k*a);
    total = total + term;
    n = n + 1;
    oldratio = newratio;
    newratio = abs(max(abs(term(:))))/abs(max(abs(total(:))));
end
TF = 1j/(k*a)^2 * total;
