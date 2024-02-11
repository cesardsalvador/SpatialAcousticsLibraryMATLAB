%% Example of Ynm
% This script plots the spherical harmonic function
% Required functions: ynm, pnm 

% Cesar D. Salvador
% daniel@si-in.com
% https://cesardsalvador.github.io/
% October 14, 2020

%% Clean and prepare the Matlab environment
close all,                                  % Close all figures
clearvars,                                  % Clear all variables in the workspace
clc,                                        % Clean the command window
workspace,                                  % Make sure the workspace panel is showing.
format longg
format compact
if(~isdeployed)
  cd(fileparts(which(mfilename)));          % Change the current folder to the folder of this m-file.
end

%% Format for vizualisation
% cmapSurf = othercolor('BuGr_14', 256);
% cmapPlot = flipud(othercolor('Greys8', 11));
cmapSurf = flipud(gray(256));
cmapPlot = gray(11);
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

%% load points from MAT files

dataFolder = 'SphericalGrid\';
subdivFactorIco = 8;
dataFileName = ['icolr', num2str(subdivFactorIco), '.mat'];
load([dataFolder, dataFileName]);
Nx = length(x);

% dataFolder = 'C:\Prog\Matlab\Data\sphgrid\md\';
% order = 5;
% Nx = (order+1)^2;
% dataFileName = ['md', sprintf('%03d', order), '.', sprintf('%05d', Nx)];
% md = load([dataFolder, dataFileName]);
% x = md(:, 1:3);
% msh = convhulln(x);
% weight = md(:, 4)/(4*pi);


%% spherical Fourier orthonormal analysis
val = 'real';
nrm = 'norm';
N = floor(sqrt(Nx))-1;

% Spherical Fourier transform matrix (Y)
Y = zeros(Nx, (N+1)^2);
for n = 0:N
    for m = -n:n
        Y(:, n^2+n+m+1) = ynm(n, m, x, val, nrm);
    end
end

% Gramm matrix, G = <Y*, Y>, to test orthonormality
G = 4*pi*Y'*diag(weight)*Y;
g = diag(G);
orthonormalityError = zeros(N+1, 1);
for n = 1:N
    indsum = n^2+1 : n^2+2*n+1;             % indices to m : |m|<=n, for every n
    num = sum(abs(g(indsum)).^2);
    den = (sum(abs(G(:, indsum)).^2, 'all') - num) / ((N+1)^2-1);
    orthonormalityError(n+1) = num/den;
end

%% Visualization
% Points on the spherical surface
figure
plot3(x(:, 1), x(:, 2), x(:, 3), '*')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
grid off
box on
axis equal
axis vis3d
titleStr = ['Spherical grid in file ', dataFileName];
title(titleStr)


% Spherical harmonic functions (Y)
figure
[azim, elev, ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
Nmax = min(2, N);
for n = 0:Nmax
    for m = -n:n
        subplot(Nmax+1, 2*Nmax+1, (2*Nmax+1)*n+Nmax+m+1)
        [v(:, 1), v(:, 2), v(:, 3)] = sph2cart(azim, elev, abs(real(Y(:, n^2+n+m+1))));
        trisurf(msh, v(:, 1), v(:, 2), v(:, 3), real(Y(:, n^2+n+m+1)), 'edgecolor', 'none')
        xlabel('$x$')
        ylabel('$y$')
        zlabel('$z$')
        grid off
        axis tight
        axis equal
        axis vis3d
        view(viewangle)
        box on
       if n == 0
           titleStr = [val, ', ', nrm, ' $Y_{nm}(\theta, \phi)$'];
           title(titleStr)
       end
    end
end

% Gramm matrix (orthonormality test)
figure
n = 0:N;
domain = [0 (N+1)^2-1];
value = log10(abs( G/max(abs(G(:))) ));
imagesc(domain, domain, value)
titleStr = [val, ', ', nrm, ' orthonormality matrix'];
title(titleStr)
xlabel('$n$')
ylabel('$n^\prime$')
set(gca, 'xtick', n.^2+n, 'xticklabel', n)
set(gca, 'ytick', n.^2+n, 'yticklabel', n)
axis xy
axis equal
axis tight
clim = [-6 0];
ctick = clim(1):clim(2);
cticklabel = ctick;
set(gca, 'clim', clim)
hcb = colorbar('eastoutside', 'ytick', ctick, 'yticklabel', cticklabel, 'ylim', clim);
clabelstr = ['$\log_{10}{| \mathbf{G} |}$, ', '$\mathbf{G}= \langle Y_{nm}, Y_{n^\prime m^\prime} \rangle$'];
set(get(hcb, 'ylabel'), 'String', clabelstr, 'interpreter', 'latex');
box on


% SNR
figure
domain = n';
value = 10*log10(abs(orthonormalityError));
plot(domain, value, '--o')
titleStr = [val, ', ', nrm, ' orthonormality error'];
title(titleStr)
xlabel('$n$')
ylabelstr = ['$10 \log{ \frac{\sum_{|m|<n}|\textrm{diag}(\mathbf{G})|^2} {\sum_{|m|<n}|\mathbf{G}-\textrm{diag}(\mathbf{G})|^2} }$'];
ylabel(ylabelstr)
set(gca, 'xtick', n, 'xticklabel', n)
xlim([domain(1) domain(end)])
ylim([min(value) max(value)])
box on

% Diagonals of Gramm matrices G1 and G2
figure
domain = 0:(N+1)^2-1;
value = abs(g);
plot(domain, value, '--o')
xlabel('$n$')
ylabelstr = ['$| \textrm{diag}(\mathbf{G}) |$'];
ylabel(ylabelstr)
set(gca, 'xtick', n.^2+n, 'xticklabel', n)
xlim([domain(1) domain(end)])
ylim([min(value) max(value)])
box on

