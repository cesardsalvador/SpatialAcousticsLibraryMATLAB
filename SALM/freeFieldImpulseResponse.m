% Set of L impulse responses of length Ns characterizing the free-field
% acoustic paths from a source position (xs) to a set of L loudspeaker positions (xl)
%
% IR = freeFieldImpulseResponse(xs, xl, Ns, Fs, c, model)
%
% INPUT:
%  xs     : Source position xs = [xs(1) xs(2) xs(3)]
%  xl     : Loudspeaker positions xl = [xl(1:L, 1) xl(1:L, 2) xl(1:L, 3)]
%  Ns     : Number of samples in time
%  Fs     : Sampling frequency
%  c      : Speed of sound in air
%  model  : Character string indicating the acoustic propagation model
%           'plane' ....................... Plane-wave model
%           'planeScaled' ................. Scaled plane-wave model used in wave field synthesis
%           'planeBandLim' ................ Band-limited plane-wave model from its truncated series expansion
%           'spherical' ................... Spherical-wave model
%           'sphericalScaled' ............. Scaled spherical-wave model used in wave field synthesis
%           'sphericalBandLim' ............ Band-limited spherical-wave model from its truncated series expansion
%           'multipoleDirection' .......... Multipole-expansion model using only directions
%           'multipoleDistance' ........... Multipole-expansion model using directions and distances
%           'virtualRecordingPlane' ......  Plane-wave virtual recording and inverse propagation
%           'virtualRecordingSpherical' ... Spherical-wave virtual recording and inverse propagation
%
% OUTPUT:
%  IR     : Impulse responses organized in a matrix of size [Ns L]
%
% See also legendrenm, besselhsph, besseljsph

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% References
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
% microphone and loudspeaker arrays”, IEEE/ACM Trans. Audio, Speech, Language Process.,
% vol. 26, no. 3, pp. 461--474, Mar. 2018.


function IR = freeFieldImpulseResponse(xs, xl, Ns, Fs, c, model)
f = (0:Ns/2)'*Fs/Ns;            % Frequency
k = 2*pi*f/c;                   % Wave number
L = size(xl, 1);                % Number of loudspeakers
N = floor(sqrt(L)-1);           % Maximum number of terms in spherical expansions
rs = rssq(xs');                 % Radial distance of source
rl = rssq(xl(1, :)');           % Radial distance of loudspeakers
cosTheta = xs*xl'/(rs*rl);      % Cosine of angle between source and loudspeakers
cosTheta(abs(cosTheta)>=1) = sign(cosTheta(abs(cosTheta)>=1));
win = cosTheta > 0;             % Spherical cap window
switch model
    case 'plane'
        TF = exp(-1j*rl*(k*cosTheta));
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'planeScaled'
        TF = exp(-1j*rl*(k*cosTheta));
        TF = -2 * TF * diag(cosTheta);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'planeBandLim'
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            jnrl = besseljsph(n, k(2:end)*rl);
            Pn = legendrenm(n, 0, cosTheta);
            termn = (1j)^n * (2*n+1) * jnrl * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        TF = conj(TF);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'spherical'
        rsl = rssq(xl-xs, 2)';
        TF = zeros(Ns/2+1, L);
        TF(2:end, :) = exp(-1j*k(2:end)*rsl) ./ (ones(size(k(2:end)))*rsl);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'sphericalScaled'
        rsl = rssq(xl-xs, 2)';
        TF = zeros(Ns/2+1, L);
        TF(2:end, :) = exp(-1j*k(2:end)*rsl);
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'sphericalBandLim'
        kind = 2;
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            hnrs = besselhsph(n, kind, k(2:end)*rs);
            jnrl = besseljsph(n, k(2:end)*rl);
            Pn = legendrenm(n, 0, cosTheta);
            termn = (-1j) * (2*n+1) * (k(2:end).*hnrs.*jnrl) * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))])) * diag(win);
    case 'multipoleDirection'
        TF = zeros(1, L);
        for n = 0:N
            Pn = legendrenm(n, 0, cosTheta);
            termn = (2*n+1) * Pn;
            TF = TF + termn;
        end
        IR = [TF; zeros(Ns-1, L)];
    case 'multipoleDistance'
        if rs < rl
            kind = 1;
        else
            kind = 2;
        end
        TF = zeros(Ns/2+1, L);
        for n = 0:N
            Pn = legendrenm(n, 0, cosTheta);
            hnrs = besselhsph(n, kind, k(2:end)*rs);
            hnrl = besselhsph(n, kind, k(2:end)*rl);
            Rnsl = hnrs./hnrl;
            regParam = rs/rl;
            Rnsl = 1./(1+regParam^2*abs(Rnsl).^2) .* Rnsl;
            termn = (2*n+1) * Rnsl * Pn;
            TF(2:end, :) = TF(2:end, :) + termn;
        end
        IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));
    otherwise
        stringDisplay = ['IR = zeros(Ns, L) because the model ', model, ...
            ' is not available. Read the HELP of freeFieldImpulseResponse.'];
        disp(stringDisplay)
        IR = zeros(Ns, L);
end
