%% Example of the spherical Fourier transform (SFT)
% Direct and inverse SFT of noisy functions on the unit sphere
% Required functions: sft, isft, ynm, pnm, randsph

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

%% Random points on the unit sphere (x)
Nx = 1000;
x = randsph(Nx, 1, 'surface');
msh = convhulln(x);
weight = quadweight(x)/(4*pi);
maxOrder = floor(sqrt(Nx))-1;

%% Noisy function on the unit sphere (s)
n = 5;
m = 2;
snrDB = 40;
snr = 10^(snrDB/10);
signal = ynm(n, m, x);
signal = signal / max(abs(signal(:)));
noise = sqrt(mean(abs(signal).^2)/snr) * randn(Nx, 1);
s = signal + noise;

%% Spherical Fourier coefficients (s1nm, s2nm and s3nm) and power densities (s1n, s2n, s3n)
s1nm = sft(s.', maxOrder, x).';
s2nm = sft(s.', maxOrder, x, 'real', 'norm', 'reg', 1/sqrt(snr)).';
s3nm = sft(s.', maxOrder, x, 'real', 'norm', 'reg', 1/sqrt(snr), weight).';
s1n = zeros(maxOrder+1, 1);
s2n = zeros(maxOrder+1, 1);
s3n = zeros(maxOrder+1, 1);
for n = 0:maxOrder
    s1n(n+1) = sum(abs(s1nm(n^2+(1:2*n+1))).^2) / (2*n+1);
    s2n(n+1) = sum(abs(s2nm(n^2+(1:2*n+1))).^2) / (2*n+1);
    s3n(n+1) = sum(abs(s3nm(n^2+(1:2*n+1))).^2) / (2*n+1);
end

%% Reconstructed signals (s1r, s2r, and s3r)
s1r = isft(s1nm.', x).';
s2r = isft(s2nm.', x).';
s3r = isft(s3nm.', x).';

%% Reconstruction errors
reconstructionErrorDB1 = 20*log10(rms(s(:)-s1r(:))/rms(s(:)))
reconstructionErrorDB2 = 20*log10(rms(s(:)-s2r(:))/rms(s(:)))
reconstructionErrorDB3 = 20*log10(rms(s(:)-s3r(:))/rms(s(:)))

%% Denoising level
denoisingDB1 = 20*log10(rms(signal(:)-s1r(:))/rms(signal(:)))
denoisingDB2 = 20*log10(rms(signal(:)-s2r(:))/rms(signal(:)))
denoisingDB3 = 20*log10(rms(signal(:)-s3r(:))/rms(signal(:)))

%% Visualization

% Points on the unit sphere
figure
plot3(x(:, 1), x(:, 2), x(:, 3), '*')
titleStr = 'Spherical random grid';
title(titleStr)
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
grid off
axis equal
axis vis3d
box on
view(viewangle)

% Functions on the unit sphere
figure
[azim, elev, ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
[v1(:, 1), v1(:, 2), v1(:, 3)] = sph2cart(azim, elev, abs(real(s)));
trisurf(msh, v1(:, 1), v1(:, 2), v1(:, 3), real(s), 'edgecolor', 'none')
titleStr = ['Polar plot of signal $s(\theta, \phi)$, SNR=', num2str(snrDB), ' dB'];
title(titleStr)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid off
axis tight
axis equal
axis vis3d
clim = [-max(abs(real(s(:)))) max(abs(real(s(:))))];
ctick = [clim(1) 0 clim(2)];
cticklabel = ctick;
set(gca, 'clim', clim)
hcb = colorbar('southoutside', 'xtick', ctick, 'xticklabel', cticklabel, 'xlim', clim);
clabelstr = '$s(\theta, \phi)$';
set(get(hcb, 'xlabel'), 'String', clabelstr, 'interpreter', 'latex');
box on
view(viewangle)

% Spherical power densitty
figure
domain = (0:maxOrder)';
value1 = 10*log10(abs(s1n));
plot(domain, value1, '*:')
hold on
value2 = 10*log10(abs(s2n));
plot(domain, value2, 'd:')
hold on
value3 = 10*log10(abs(s3n));
plot(domain, value3, 'o:')
titleStr = ['Power density of $s_{nm}=\textrm{SFT}\{s(\theta, \phi)\}$, SNR=', num2str(snrDB), ' dB'];
title(titleStr)
xlabel('$n$')
ylabel('$10 \log_{10}\left|\frac{1}{2n+1}\sum_{|m| \leq n}|s_{nm}|^2\right|$')
xlim([min(domain) max(domain)])
ylim([min([value1(:); value2(:); value3(:)]) max([value1(:); value2(:); value3(:)])])
grid on
box on
legend('Moore-Penrose', 'Tikhonov', 'Tikhonov, quad', 'location', 'northeast')

% Reconstructed functions on the unit sphere
figure
subplot(131)
[azim, elev, ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
[v1(:, 1), v1(:, 2), v1(:, 3)] = sph2cart(azim, elev, abs(real(s1r)));
trisurf(msh, v1(:, 1), v1(:, 2), v1(:, 3), real(s1r), 'edgecolor', 'none')
titleStr = 'Moore-Penrose';
title(titleStr)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid off
axis tight
axis equal
axis vis3d
clim = [-max(abs(real(s1r(:)))) max(abs(real(s1r(:))))];
ctick = [clim(1) 0 clim(2)];
cticklabel = ctick;
set(gca, 'clim', clim)
hcb = colorbar('southoutside', 'xtick', ctick, 'xticklabel', cticklabel, 'xlim', clim);
clabelstr = '$\hat{s}(\theta, \phi)$';
set(get(hcb, 'xlabel'), 'String', clabelstr, 'interpreter', 'latex');
box on
view(viewangle)

subplot(132)
[azim, elev, ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
[v2(:, 1), v2(:, 2), v2(:, 3)] = sph2cart(azim, elev, abs(real(s2r)));
trisurf(msh, v2(:, 1), v2(:, 2), v2(:, 3), real(s2r), 'edgecolor', 'none')
titleStr = 'Tikhonov';
title(titleStr)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid off
axis tight
axis equal
axis vis3d
clim = [-max(abs(real(s2r(:)))) max(abs(real(s2r(:))))];
ctick = [clim(1) 0 clim(2)];
cticklabel = ctick;
set(gca, 'clim', clim)
hcb = colorbar('southoutside', 'xtick', ctick, 'xticklabel', cticklabel, 'xlim', clim);
clabelstr = '$\hat{s}(\theta, \phi)$';
set(get(hcb, 'xlabel'), 'String', clabelstr, 'interpreter', 'latex');
box on
view(viewangle)

subplot(133)
[azim, elev, ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
[v2(:, 1), v2(:, 2), v2(:, 3)] = sph2cart(azim, elev, abs(real(s3r)));
trisurf(msh, v2(:, 1), v2(:, 2), v2(:, 3), real(s3r), 'edgecolor', 'none')
titleStr = 'Tikhonov, quad';
title(titleStr)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid off
axis tight
axis equal
axis vis3d
clim = [-max(abs(real(s3r(:)))) max(abs(real(s3r(:))))];
ctick = [clim(1) 0 clim(2)];
cticklabel = ctick;
set(gca, 'clim', clim)
hcb = colorbar('southoutside', 'xtick', ctick, 'xticklabel', cticklabel, 'xlim', clim);
clabelstr = '$\hat{s}(\theta, \phi)$';
set(get(hcb, 'xlabel'), 'String', clabelstr, 'interpreter', 'latex');
box on
view(viewangle)

