% Total pressure over a rigid sphere (total = incident + scattered)
%   p = pressureRigidSphere(a, b, theta, k, thres)
%
% INPUT
% a:        radius of the point of measurement on the rigid spherical surface
% b:        radius of the radiating point source
% theta:    incident angle between the point of measurement and the point source in DEGREES
% k:        wave number (k = 2*pi*frequency/velocity_of_sound_in_air)
% thres:    threshold error (typical 0.0001)
%
% OUPUT
% p:        total pressure on the point of measurement
%
% See also sphere_transfer_function.

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


function p = pressureRigidSphere(a, b, theta, k, thres)

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

p = (rho * exp(-1j * mu) * total) ./ (1j * mu);
p = p .* exp(1j*k*b)/b;
