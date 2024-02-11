% This function creates a template for a Matlab script in the current folder
%
% >> newScript(fileName)
%
% INPUT:
%  fileName: Character string containing the file name (eg. 'myScript.m')
%
% See also newFunction

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [] = newScript(fileName)
folderNameTemplate = 'C:\Prog\Matlab\SpatialAcousticLibrary\';              % Set this folder to a local path
fileNameTemplate = 'templateScript.m';                       % Verify that this template is in the local path
copyfile([folderNameTemplate, fileNameTemplate], fileName)   % Creates a copy in the current folder
edit(fileName)                                               % Opens the created copy
end
