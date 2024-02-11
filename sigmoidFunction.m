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

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% June 17, 2021

% REFERENCES
% [1] E. G. Williams, Fourier Acoustics: Sound Radiation and Nearfield Acoustical Holography,
% London, UK: Academic Press, 1999.
% [2] C. D. Salvador et al., “Boundary matching filters for spherical microphone and loudspeaker arrays”,
% IEEE/ACM Trans. Audio, Speech, Language Process., vol. 26, no. 3, pp. 461--474, Mar. 2018,
% [3] Silicon Integrated Co., Ltd., https://www.si-in.com/

function y = sigmoidFunction(x, scale, shift)
y = 1./(1 + exp(-scale*(x-shift)));
