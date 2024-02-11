% Example of sphere point picking using the function RANDSPH

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% October 14, 2020

close all
clearvars
clc

%% Parameters for visualization
cmapSurf = parula(256);
cmapPlot = gray(256);
font_size = 18;
set(0, 'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmapSurf, ...
    'defaultAxesColorOrder', cmapPlot, ...
    'defaultAxesFontName', 'times', ...
    'defaultTextFontSize', font_size, ...
    'defaultAxesFontSize', font_size, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')

%% Parameters for figure files
figFolder = 'Fig/';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';
%figFormat = '-dmeta'; figNameExt = '.emf'; figRender = 'zbuffer'; figRes = '-r150';

%% load points from MAT files
dataFolder = 'SphericalGrid\';
subdivFactorIco = 14;
dataFileName = ['icolr', num2str(subdivFactorIco), '.mat'];
load([dataFolder, dataFileName]);
Nx = length(x);
sph.xico = x;
sph.xrand = randsph(Nx, 1, 'surface');

%% spherical Fourier orthonormal analysis
val = 'complex';
nrm = 'sch';
maxOrder = floor(sqrt(Nx))-1;
Yico = zeros(Nx, (maxOrder+1)^2);
for n = 0:maxOrder
    for m = -n:n
        Yico(:, n^2+n+m+1) = ynm(n, m, x, val, nrm);
    end
end

%% Visualization
figure
plot3(sph.xico(:, 1), sph.xico(:, 2), sph.xico(:, 3), 'b.')
hold on
plot3(sph.xrand(:, 1), sph.xrand(:, 2), sph.xrand(:, 3), 'r.')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
grid off
box on
axis equal
axis vis3d

