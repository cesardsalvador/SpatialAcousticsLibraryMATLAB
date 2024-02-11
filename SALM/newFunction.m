% This function creates a template for a Matlab function in the current folder
%
% >> newFunction(fileName)
%
% INPUT:
%  fileName: Character string containing the file name (eg. 'myScript.m')
%
% See also newScript

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [] = newFunction(fileName)
folderNameTemplate = 'C:\Prog\Matlab\SpatialAcousticLibrary\';              % Set this folder to a local path
fileNameTemplate = 'templateFunction.m';                     % Verify that this template is in the local path
copyfile([folderNameTemplate, fileNameTemplate], fileName)   % Creates a copy in the current folder
edit(fileName)                                               % Opens the created copy
end
