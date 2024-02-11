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

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% Reference and citation
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.
%     DOI: 10.1109/TASLP.2017.2778562
% [2] C. D. Salvador et al., “Design theory for binaural synthesis:
%     Combining microphone array recordings and head-related transfer function datasets,”
%     Acoust. Sci. Technol., vol. 38, no. 2, pp. 51–62, Mar. 2017.
%     DOI: 10.1250/ast.38.51


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
