% Impulse response of a feedback delay network with L delay lines
%
%  IR = feedbackDelayImpulseResponse(delay, feedbackMatrix, inputGain, directGain, outputGain, Ns, Fs)
%
% INPUT:
%  inputGain            : Column vector with L input gains
%  outputGain           : Column vector with L output gains
%  directGain           : Scalar gain of the direct path
%  delay                : Column vector with L delay values in samples
%  feedbackMatrix       : Square matrix of size L with the scalar feedback coefficients
%  Ns                   : Number of samples of the impulse response
%  Fs                   : Sampling frequency in Hertz
%
% OUTPUT:
%  IR                   : Column vector containing the impulse response of Ns samples
%
% See also fft, ifft

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% REFERENCES
% [1] S. J. Schlecht et al., “Scattering in feedback delay networks,”
% IEEE/ACM Trans. Audio, Speech, Language Process., vol. 28, pp. 1915–1924, 2020.

function IR = feedbackDelayImpulseResponse(delay, feedbackMatrix, inputGain, directGain, outputGain, Ns, Fs)
f = (0:Ns/2)'*Fs/Ns;                                    % Frequency domain
z = exp(1j*2*pi*f/Fs);                                  % Unit circle on the z-domain
TF = zeros(Ns/2+1, 1);
for I = 1:Ns/2+1
    P = diag(z(I).^delay) - feedbackMatrix;             % Polynomial matrix or loop transfer function
    TF(I) = outputGain'*(P\inputGain) + directGain;     % Transfer function of the network
end
IR = real(ifft([TF; conj(TF(end-1:-1:2))]));            % Impulse response of the network
