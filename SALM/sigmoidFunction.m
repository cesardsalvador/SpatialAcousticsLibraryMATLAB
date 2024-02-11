% This functions applies a gain to an input signal and outputs the result
%
%  y = sigmoidFunction(x, scale, shift)
%
% INPUT:
%  x            : Domain
%  scale        : Scale to control the slope and inflection
%
% OUTPUT:
%  y            : Sigmoid function of x
%
% See also exampleSigmoidFunction.m

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


function y = sigmoidFunction(x, scale, shift)
y = 1./(1 + exp(-scale*(x-shift)));
