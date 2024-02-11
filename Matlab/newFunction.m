% This function creates a template for a Matlab function in the current folder
%
% >> newFunction(fileName)
%
% INPUT:
%  fileName: Character string containing the file name (eg. 'myScript.m')
%
% See also newScript

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 6, 2020

% REFERENCE
% https://www.mathworks.com/matlabcentral/answers/56672-how-can-i-create-a-standard-matlab-template-for-new-programs

function [] = newFunction(fileName)
folderNameTemplate = 'C:\Prog\Matlab\SpatialAcousticLibrary\';              % Set this folder to a local path
fileNameTemplate = 'templateFunction.m';                     % Verify that this template is in the local path
copyfile([folderNameTemplate, fileNameTemplate], fileName)   % Creates a copy in the current folder
edit(fileName)                                               % Opens the created copy
end
