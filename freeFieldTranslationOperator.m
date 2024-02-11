% Free-field translation operator for sound pressure fields.
%
%  [IR, TF] = freeFieldtranslationOperator(xSource, xInitial, xFinal, k, model, direction)
%
% INPUT:
%   xSource     : Position of sound source in Cartesian coordinates ([x1 y1 z1; ...; xP yP zP]).
%   xInitial	: Position of measurement in Cartesian coordinates ([x y z]).
%   xFinal      : Position of translation in Cartesian coordinates ([x y z]).
%   k           : Wave number. Column vector of length Ns/2+1, where Ns is number of samples in time.
%   model       : Sound propagation model of translation operator:
%                 'plane'       :Plane-wave propagation model.
%                 'spherical'   :Spherical-wave propagation model.
%   direction   : Direction of translation:
%                 'direct'      : Direct translation from xInitial to xFinal
%                 'inverse'     : Inverse translation from xFinal to xInitial
%
% OUTPUT:
%   IR          : Inpulse response of translation operator. Matrix of size [Ns P].
%   TF          : Transfer function of translation operator. . Matrix of size [Ns/2+1 P].
%
% See also exampleFreeFieldTranslationOperator.m, exampleDistanceVaryingFilterSphericalEarCentering.m

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% January 4, 2021

% REFERENCES
% [1] E. G. Williams, Fourier Acoustics: Sound Radiation and Nearfield Acoustical Holography,
% London, UK: Academic Press, 1999.
% [2] I. Ben Hagai et al., “Acoustic centering of sources measured by surrounding
% spherical microphone arrays,” J. Acoust. Soc. Am., vol. 130, no. 4, pp. 2003–2015, 2011.
% [3] N. R. Shabtai et al., “Acoustic centering of sources with high-order radiation patterns,”
% J. Acoust. Soc. Am., vol. 137, no. 4, pp. 1947–1961, 2015.
% [4] Z. Ben-Hur et al., “Efficient representation and sparse sampling of head-related transfer functions
% using phase-correction based on ear alignment,” IEEE Trans. Audio, Speech, Language Process., pp. 1–1, 2019.

function [IR, TF] = freeFieldTranslationOperator(xSource, xInitial, xFinal, k, model, direction)
Ns = 2*(length(k)-1);                                   % Number of samples in time.
xI2S = xSource - xInitial;                              % Relative position of source w.r.t. initial position.
distI2S = rssq(xI2S, 2);                                % Distance from source to initial position.
switch model
    case 'plane'
        xI2F = xFinal - xInitial;                       % Relative position of final w.r.t. initial position.
        distI2F = rssq(xI2F, 2);                        % Distance from final to initial position.
        cosAngleSIF = xI2S*xI2F'./(distI2S*distI2F);    % Cosine of angle between source and final relative positions.
        directTF = exp(-1j*k*distI2F*cosAngleSIF');     % Direct plane-wave transfer function of translator.
        if strcmp(direction, 'direct')
            TF = directTF;
        elseif strcmp(direction, 'inverse')
            TF = 1./directTF;                           % Inverse plane-wave transfer function of translator.
        else
            TF = zeros(Ns/2+1, 1);
            disp('Output zero. DIRECTION must be DIRECT or INVERSE.')
        end
    case 'spherical'
        distF2S = rssq(xSource - xFinal, 2);            % Distance from source to final position.
        directTF = exp(-1j*k*(distF2S-distI2S)') * ...
                   diag(distI2S./distF2S);              % Direct spherical-wave transfer function.
        if strcmp(direction, 'direct')
            TF = directTF;
        elseif strcmp(direction, 'inverse')
            TF = 1./directTF;                           % Inverse spherical-wave transfer function of translator.
        else
            TF = zeros(Ns/2+1, 1);
            disp('Output zero. DIRECTION must be DIRECT or INVERSE.')
        end
    otherwise
        TF = zeros(Ns/2+1, 1);
        disp('Output zero. TYPE must be PLANE or SPHERICAL.')
end
IR = real(ifft([TF; conj(TF(end-1:-1:2, :))]));         % Impulse response of translator.
