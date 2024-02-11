% This functions applies a gain to an input signal and outputs the result
%
%  outputSignal = templateFunction(inputSignal, gain)
%
% INPUT:
%  inputSignal	: Input signal
%  gain         : Gain
%
% OUTPUT:
%  outputSignal	: Output signal
%
% See also FUNCTION1, FUNCTION2, ...

% Cesar D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% March 5, 2023

% REFERENCES
% [1] E. G. Williams, Fourier Acoustics: Sound Radiation and Nearfield Acoustical Holography,
% London, UK: Academic Press, 1999.
% [2] C. D. Salvador et al., “Boundary matching filters for spherical microphone and loudspeaker arrays”,
% IEEE/ACM Trans. Audio, Speech, Language Process., vol. 26, no. 3, pp. 461--474, Mar. 2018,
% [3] Silicon Integrated Co., Ltd., https://www.si-in.com/

function outputSignal = templateFunction(inputSignal, gain)
outputSignal = gain * inputSignal;
