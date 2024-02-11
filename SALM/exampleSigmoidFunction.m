%% Example of sigmoid function

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
cmapSurf = jet(256);
cmapLine = flipud(cool(2));
fontSize = 16;
lineWidth = 4;
set(0, ...
    'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmapSurf, ...
    'defaultAxesColorOrder', cmapLine, ...
    'DefaultLineLineWidth', lineWidth, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', fontSize, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')

%% Format for saving figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Main processing
Ns = 1024;                                  % Number of samples
Fs = 16e3;                                  % Sampling frequency
t = (0:Ns-1)' * 1/Fs;                       % Time domain
shift = 30e-3;                             % Transition center time in seconds
transitionDuration = 25e-3;                % Transition duration in seconds
scale = 10/transitionDuration;
windowLate = sigmoidFunction(t, scale, shift);
windowEarly = sigmoidFunction(t, -scale, shift);

%% Visualization
figure
domain = t*1e3;                             % Time domain in ms (x axis)
value = windowEarly;                        % Early window (y axis)
plot(domain, value)
hold on
value = windowLate;                         % Late window (y axis)
plot(domain, value)
title('Sigmoid windows')
xlabel('$t$ (ms)')
ylabel('$w(t)$')
xlim([min(domain) max(domain)])
grid on
box on

%% Saving figures
% figName = ['figFourierTransform', num2str(f0), 'HzTone'];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

