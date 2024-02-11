% This function calculates the impulse response (IR) and transfer function (TF) of a
% diffuse-field filter to equalize a dataset of raw acoustic impulse responses.
%
%  [equalizerIR, equalizerTF] = diffuseFieldFilter(datasetIR, phaseType)
%  [equalizerIR, equalizerTF] = diffuseFieldFilter(datasetIR, phaseType, weight)
%
% INPUT:
%  datasetIR    : Matrix of size [Ns L] with Ns-samples impulse responses for L positions
%  phaseType    : 'zero' or 'minimum' to indicate the type of phase of the equalizer.
%  weight       : Vector of size [L 1] with integration weights s.t. sum(weight)=1
%
% OUTPUT:
%  equalizerIR  : Vector of size [Ns 1] with the equalizer impulse response
%  equalizerTF  : Vector of size [Ns/2+1 1] with the equalizer transfer function
%
% See also exampleDiffuseFieldFilter

% December 15, 2020
% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/

% REFERENCES
% [1] J. C. Middlebrooks, “Individual differences in external-ear transfer functions reduced
% by scaling in frequency,” J. Acoust. Soc. Am., vol. 106, no. 3, pp. 1480–1492, 1999.
% [2] J. C. Middlebrooks et al., “Directional dependence of interaural envelope delays,”
% J. Acoust. Soc. Am., vol. 87, no. 5, pp. 2149–2162, 1990.

function [equalizerIR, equalizerTF] = diffuseFieldFilter(datasetIR, phaseType, varargin)
[Ns, L] = size(datasetIR);
datasetTF = fft(datasetIR);                                         % Dataset transfer functions for L positions
datasetTF = datasetTF(1:Ns/2+1, :);
if nargin < 3
    magnitude = 1 ./ sqrt(sum(abs(datasetTF).^2 / L, 2));           % Magnitude of equalizer
else
    magnitude = 1 ./ sqrt(sum(abs(datasetTF).^2 * ...
                          diag(varargin{1}), 2));                   % Magnitude of equalizer with integration weights
end
switch phaseType
    case 'zero'
        phase = zeros(Ns/2+1, 1);                                   % Zero phase
    case 'minimum'
        phase = imag(hilbert(log([1./magnitude; ...
                                  1./magnitude(end-1:-1:2)])));     % Minimum phase
        phase = phase(1:Ns/2+1);
    otherwise
        phase = zeros(Ns/2+1, 1);
        disp('PHASETYPE not available. Phase = 0.')
end
equalizerTF = magnitude .* exp(1j*phase);                           % Transfer function of equalizer
equalizerIR = real(ifft([equalizerTF; ...
                         conj(equalizerTF(end-1:-1:2))]));          % Impulse response of equalizer
