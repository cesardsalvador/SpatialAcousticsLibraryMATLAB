%% Example of freeFieldTranslationOperator
% This script plots the Fourier transform of a single-frequency tone
% This script requires the following libraries:
% LIBRARY1 available at HYPERLINK1
% LIBRARY2 available at HYPERLINK2
% ...

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


%% Clean and prepare the Matlab environment
close all                                   % Close all figures
clearvars                                   % Clear all variables in the workspace
clc                                         % Clean the command window
workspace                                   % Make sure the workspace panel is showing.
format longg
format compact

%% Format for vizualisation
%cmap = othercolor('BuDRd_18', 256);        % https://www.mathworks.com/matlabcentral/fileexchange/30564-othercolor
cmap = jet(256);
fontSize = 18;
set(0, ...
    'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmap, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', fontSize, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')
viewangle = [130 28];

%% Format for saving figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Time (t), frequency (f) and wave-number (k) domains
Ns = 512;                                   % Number of samples
Fs = 48e3;                                  % Sampling frequency
c = 343;                                    % Speed of sound in air
t = (0:Ns-1)' * 1/Fs;                       % Time (seconds)
f = (0:Ns/2)' * Fs/Ns;                      % Frequency (Hertz)
k = 2*pi*f/c;                               % Wave number (1/meter)

%% Geometry
xSource = [1 1 1];                          % Position of sound source in Cartesian coordinates.
xInitial = [0 .5 0];                        % Position of measurement in Cartesian coordinates.
xFinal = [0 0 0];                           % Position of translation in Cartesian coordinates.

%% Free-field translation opetator: Impulse response (IR) and transfer function (TF)
model = 'spherical';                        % Sound propagation model of translation operator.
direction = 'direct';                       % Direction of translation.
[IR, ~] = freeFieldTranslationOperator(xSource, xInitial, xFinal, k, model, direction);

%% Visualization
figure
domain = t*1e3;                             % Time domain in miliseconds (x axis)
value = IR;                                 % Value of signal (y axis)
plot(domain, value, 'b', ...
    'lineWidth', 2)
title('Impulse response of translation operator')
xlabel('$t$ (ms)')
ylabel('$h(t)$')
xlim([min(domain) max(domain)])
ylim([min(value) max(value)])
grid on
box on

%% Saving figures
% figName = ['figFourierTransform', num2str(f0), 'HzTone'];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

