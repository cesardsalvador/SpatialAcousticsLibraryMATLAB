% Set of L impulse responses of length Ns characterizing the free-field
% acoustic paths from a source position (xs) to a set of L loudspeaker positions (xl)
%
% IR = freeFieldImpulseResponse(xs, xl, k, model)
%
% INPUT:
%  xs     : Source position xs = [xs(1) xs(2) xs(3)]
%  xl     : Loudspeaker positions xl = [xl(1:L, 1) xl(1:L, 2) xl(1:L, 3)]
%  k      : Wave number k = [k(1); k(2); ... ; k(Ns/2+1)], where Ns is the frame length
%  model  : Character string indicating the acoustic propagation model
%           'plane' ....................... Plane-wave model
%           'planeBandLim' ................ Plane-wave model with limited spatial bandwidth
%           'spherical' ................... Spherical-wave model
%           'sphericalBandLim' ............ Spherical-wave model with limited spatial bandwidth
%           'multipoleDirection' .......... Multipole-expansion model using only directions
%           'multipoleDistance' ........... Multipole-expansion model using directions and distances
%           'virtualRecordingPlane' ......  Plane-wave virtual recording and inverse propagation
%           'virtualRecordingSpherical' ... Spherical-wave virtual recording and inverse propagation
%
% OUTPUT:
%  IR     : Impulse responses organized in a matrix of size [Ns L]
%
% See also pnm, besselhsph, besseljsph

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% References
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
% microphone and loudspeaker arrays”, IEEE/ACM Trans. Audio, Speech, Language Process.,
% vol. 26, no. 3, pp. 461--474, Mar. 2018.

function IR = freeFieldImpulseResponseBackUp(xs, xl, k, model)
Ns = 2*(length(k)-1);           % Number of samples in time
L = size(xl, 1);                % Number of loudspeakers
N = floor(sqrt(L)-1);           % Maximum number of terms in spherical expansions
rs = rssq(xs');                 % Radial distance of source
rl = rssq(xl(1, :)');           % Radial distance of loudspeakers
cosTheta = xs*xl'/(rs*rl);      % Cosine of angle between source and loudspeakers
cosTheta(abs(cosTheta)>=1) = sign(cosTheta(abs(cosTheta)>=1));
% win = cosTheta > 0.1;           % spherical cap window
switch model
    case 'plane'
        TF = exp(-1j*rl*(k*cosTheta));
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    case 'spherical'
        rsl = rssq(xl-xs, 2)';
        TF = zeros(Ns/2+1, L);
        TF(2:end, :) = exp(-1j*k(2:end)*rsl) ./ (ones(size(k(2:end)))*rsl);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    case 'planeBandLim'
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            jnrl = besseljsph(n, k(2:end)*rl);
            Pn = pnm(n, 0, cosTheta);
            termn = (1j)^n * (2*n+1) * jnrl * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        TF = conj(TF);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    case 'sphericalBandLim'
        kind = 2;
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            hnrs = besselhsph(n, kind, k(2:end)*rs);
            jnrl = besseljsph(n, k(2:end)*rl);
            Pn = pnm(n, 0, cosTheta);
            termn = (-1j) * (2*n+1) * (k(2:end).*hnrs.*jnrl) * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    case 'multipoleDirection'
        TF = zeros(1, L);
        for n = 0:N
            Pn = pnm(n, 0, cosTheta);
            termn = (2*n+1) * Pn;
            TF = TF + termn;
        end
        IR = [TF; zeros(Ns-1, L)];
    case 'multipoleDistance'
        kind = 1;
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            hnrs = besselhsph(n, kind, k(2:end)*rs) ./ besselhsph(0, kind, k(2:end)*rs);
            hnrl = besselhsph(n, kind, k(2:end)*rl) ./ besselhsph(0, kind, k(2:end)*rl);
            Pn = pnm(n, 0, cosTheta);
            termn = (2*n+1) * (hnrl./hnrs) * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    otherwise
        stringDisplay = ['IR = zeros(Ns, L) because the model ', model,' is not available. Read the HELP of freeFieldImpulseResponse.'];
        disp(stringDisplay)
        IR = zeros(Ns, L);
end
