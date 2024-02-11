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

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function outputSignal = templateFunction(inputSignal, gain)
outputSignal = gain * inputSignal;
