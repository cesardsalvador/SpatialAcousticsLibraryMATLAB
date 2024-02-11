%% Example of Pnm
% This script plots the associated Legendre polynomial
% Required function: pnm 

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


%% Clean and prepare the Matlab environment
close all,                                  % Close all figures
clearvars,                                  % Clear all variables in the workspace
clc,                                        % Clean the command window
workspace,                                  % Make sure the workspace panel is showing.
format longg
format compact

%% Format for vizualisation
cmapSurf = flipud(gray(256));
cmapPlot = gray(10);
font_size = 14;
set(0, 'defaultFigureColor', 'w', ...
    'defaultLineLineWidth', 1, ...
    'defaultFigureColorMap', cmapSurf, ...
    'defaultAxesColorOrder', cmapPlot, ...
    'defaultAxesFontName', 'times', ...
    'defaultTextFontSize', font_size, ...
    'defaultAxesFontSize', font_size, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')
viewangle = [130 28];

%% Format for saving figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Associated Legendre polinomial
nrm = 'norm';
N = 14;
x = -1:.01:1;
Nx = length(x);
p1 = zeros(Nx, N+1, N+1);
p2 = p1;
for n = 0:N
    p = legendre(n, x, nrm);
    for m = 0:n
        p1(:, n+1, m+1) = p(m+1, :);
        p2(:, n+1, m+1) = pnm(n, m, x, nrm);
    end
end
p1 = reshape(p1, [Nx (N+1)^2]);
p2 = reshape(p2, [Nx (N+1)^2]);
error = max(abs(p1(:)-p2(:)))

%% Visualization
figure
domain = x;
value = 20*log10(abs(p2));
plot(domain, value, ':', 'lineWidth', .7, 'color', 'k')
title(['Associated Legendre polynomial'])
xlabel('$x$')
ylabel('$20 \log_{10}|P_n^m(x)|$')
xlim([min(domain) max(domain)])
ylim([-6 6])
set(gca, 'xscale', 'linear', 'yscale', 'linear')
%set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
%set(gca, 'ytick', -60:6:60, 'yticklabel', -60:6:60)
grid on
box on

%% Saving figures
% figName = ['figAssociatedLegendrePolynomial'];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

