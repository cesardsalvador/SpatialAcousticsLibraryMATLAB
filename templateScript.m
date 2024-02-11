%% Description of script
% This script plots the Fourier transform of a single-frequency tone
% This script requires the following libraries:
% LIBRARY1 available at HYPERLINK1
% LIBRARY2 available at HYPERLINK2
% ...

% Cesar D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% March 5, 2023

%% Clean and prepare the Matlab environment
close all                                   % Close all figures
clearvars                                   % Clear all variables in the workspace
clc                                         % Clean the command window
workspace                                   % Make sure the workspace panel is showing

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

%% Parameters to export figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Main processing
Ns = 512;                                   % Number of samples
Fs = 48e3;                                  % Sampling frequency
t = (0:Ns-1)' * 1/Fs;                       % Time domain
f = (0:Ns/2)' * Fs/Ns;                      % Frequency domain
f0 = 1000;                                  % Frequency of signal
s = cos(2*pi*f0*t);                         % Signal in time domain
S = fft(s); S = S(1:Ns/2+1);                % Signal in frequency domain

%% Visualization
figure
domain = f*1e-3;                            % Frequency domain in kHz (x axis)
value = 20*log10(abs(S));                   % Value of signal (y axis)
plot(domain, value, 'b', ...
    'lineWidth', 2)
title(['Fourier transform of a ', sprintf('%.2f', f0), ' Hz tone'])
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|S(f)|$')
xlim([min(domain) max(domain)])
ylim([min(value) max(value)])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
set(gca, 'ytick', -60:6:60, 'yticklabel', -60:6:60)
grid on
box on

%% Exporting figures
% figName = ['figFourierTransform', num2str(f0), 'HzTone'];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

