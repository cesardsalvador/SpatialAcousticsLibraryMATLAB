%% Example of spherical distance-varying filters (DVF)
%
% Summary:
%   1) Reads a dense spherical dataset of HRIRs from a SOFA file.
%   2) Applies diffuse field equalization to the dense HRIR dataset.
%   3) Select a sparse spherical subset of HRIRs from the dense set.
%   4) Performs distance extrapolation on the sparse spherical set
%      from a original distance "a" to a new distance "b".
%
% Notation:
%   h: head-related impulse response (HRIR)
%   H: head-related transfer function (HRTF)
%   hnm: spherical Fourier transform of h
%   Hnm: spherical Fourier transform of H
%
% Requirements:
%   1) SOFA Matlab/Octave API version 1.0.2
%      https://sourceforge.net/projects/sofacoustics/
%   2) Near-distance HRTF dataset
%      https://cesardsalvador.github.io/download.html
%
% December 16, 2020
% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
%
% REFERENCES
% [1] C. D. Salvador et al., “Dataset of near-distance head-related transfer functions
% calculated using the boundary element method,” Proc. Audio Eng. Soc. Int. Conf.
% Spatial Reproduction —Aesthetics and Science—, Tokyo, Japan, Aug. 2018.
% [2] R. Duraiswami et al., “Interpolation and range extrapolation of HRTFs,”
% Proc. IEEE ICASSP, May 2004, vol. 4, pp. 45–48.
% [3] J. C. Middlebrooks et al., “Directional dependence of interaural envelope delays,”
% J. Acoust. Soc. Am., vol. 87, no. 5, pp. 2149–2162, 1990.

%% Clean and prepare the Matlab environment
close all                                   % Close all figures
clearvars                                   % Clear all variables in the workspace
clc                                         % Clean the command window
workspace                                   % Make sure the workspace panel is showing.
format longg
format compact

%% Parameters for visualization
cmap = parula(256);
font_size = 14;
set(0, 'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmap, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', font_size, ...
    'defaultAxesFontSize', font_size, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')
viewAngle = [-40 10];

%% Parameters for figure files
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r150';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Sofa files (HRIR)
SOFAstart;
listenerModelType = 'ht';                              % 'h' (only head) or 'ht' (head and torso)
%sourceDistribution = 'sph_ico_dir642_dist93';
sourceDistribution = 'sph_ico_dir2562_dist21';         % Dense directions are useful for diffuse-field equalization
%sourceDistribution = 'sph_md_dir625_dist93';
%sourceDistribution = 'sph_md_dir2500_dist21'; integration_weight_file = 'integration_weight_md_dir2500.mat';

%% Structure array Obj from SOFA file (HRIR)
sofaFolder = 'C:\Prog\Matlab\Data\hrir\sofa\distance\';         % Set your local folder
sofaFileName = dir([sofaFolder, '*', listenerModelType, '_',sourceDistribution, '.sofa']);
indHead = 1;                                                    % 1:bkhats, 2:samrai, 3:sakamoto, 4:han
Obj = SOFAload([sofaFolder, sofaFileName(indHead).name]);
Obj = orderfields(Obj);

%% Dimensions of data (Structure array Obj)
FsFilter = Obj.Data.SamplingRate;                    % Sampling frequency in Hertz
Ns = Obj.API.N;                                      % Number of samples along time describing a measurement for one position
sampleCircshift = Obj.Data.Delay;                    % Number of samples used in circular shift of HRIRs (circular delay)
Ndir = Obj.SourcePositionNumberDirections;           % Number of source directions
Ndist = Obj.SourcePositionNumberDistances;           % Number of source distances
Npos = Obj.API.M;                                    % Number of source positions around the listener (Npos = Ndir * Ndist)
c = Obj.SpeedOfSoundInMPS;                           % Speed of sound in air (m/s)

%% HRIR dense set: hLeftDense and hRightDense (size: Ns * Ndir * Ndist)
headModelName = Obj.GLOBAL_ListenerShortName;
xEarLeft = Obj.ReceiverPosition(1, :);
xEarRight = Obj.ReceiverPosition(2, :);
hLeftDense = permute(reshape(squeeze(Obj.Data.IR(:, 1, :)), [Ndir Ndist Ns]), [3 1 2]);
hRightDense = permute(reshape(squeeze(Obj.Data.IR(:, 2, :)), [Ndir Ndist Ns]), [3 1 2]);

%% Cirular shift of HRIR dense set to fully include them within a single frame
circularShiftSamples = -165;
hLeftDense = circshift(hLeftDense, circularShiftSamples);
hRightDense = circshift(hRightDense, circularShiftSamples);

%% One-dimensional domains
t = (0:Ns-1)'*1/FsFilter;                               % Time (seconds)
f = (0:Ns/2)'*FsFilter/Ns;                              % Frequency (Hertz)
c = 344;                                                % Speed of sound in air
k = 2*pi*f/c;                                           % Wave number

%% Source positions in spherical and cartesian coordinates
azimh = Obj.SourcePosition(:, 1);
elevh = Obj.SourcePosition(:, 2);
disth = Obj.SourcePosition(:, 3);
[xh(:, 1), xh(:, 2), xh(:, 3)] = sph2cart(azimh*pi/180, elevh*pi/180, disth);
azimh = reshape(azimh, [Ndir Ndist]);                     % Azimuthal angle in degree (size: Ndir * Ndist)
elevh = reshape(elevh, [Ndir Ndist]);                     % Elevation angle in degree (size: Ndir * Ndist)
disth = reshape(disth, [Ndir Ndist]);                     % Radial distance in metre (size: Ndir * Ndist)
xh = reshape(xh, [Ndir Ndist 3]);                         % Cartesian coordinates in metre (size: Ndir * Ndist * 3)

%% Select HRIR for a dense spherical grid (xh) at a single distance
distanceCm = 200;                                 % Distance (cm) from head center to point source (200, 150, 100:-1:10)
[~, indDist] = min(abs(disth(1,:)-distanceCm*1e-2));
azimh = azimh(:, indDist);
elevh = elevh(:, indDist);
disth = disth(:, indDist);
xh = squeeze(xh(:, indDist, :));
hLeftDenseSingleDistance = hLeftDense(:, :, indDist);
hRightDenseSingleDistance = hRightDense(:, :, indDist);

%% HRTF (H): FFT of HRIR (h)
HLeftDenseSingleDistance = fft(hLeftDenseSingleDistance);
HLeftDenseSingleDistance = HLeftDenseSingleDistance(1:Ns/2+1, :);
HRightDenseSingleDistance = fft(hRightDenseSingleDistance);
HRightDenseSingleDistance = HRightDenseSingleDistance(1:Ns/2+1, :);

%% Quadrature weights for integration on the dense unit sphere
[xhUnitSphere(:, 1), xhUnitSphere(:, 2), xhUnitSphere(:, 3)] = sph2cart(azimh*pi/180, elevh*pi/180, 1);
weightDense = quadweight(xhUnitSphere)/(4*pi);

%% Diffuse-field equalizer (equalizer: impulse response, Equalizer: transfer function)
phaseType = 'zero';                 % 'zero' or 'minimum' to indicate the type of phase of the equalizer
[equalizerLeft, EqualizerLeft] = diffuseFieldFilter(hLeftDenseSingleDistance, phaseType, weightDense);
[equalizerRight, EqualizerRight] = diffuseFieldFilter(hRightDenseSingleDistance, phaseType, weightDense);

%% Diffuse-field equalization
HLeftDenseSingleDistanceEqualized = diag(EqualizerLeft) * HLeftDenseSingleDistance;
HRightDenseSingleDistanceEqualized = diag(EqualizerRight) * HRightDenseSingleDistance;
hLeftDenseSingleDistanceEqualized = real(ifft([HLeftDenseSingleDistanceEqualized; ...
    conj(HLeftDenseSingleDistanceEqualized(end-1:-1:2, :))]));
hRightDenseSingleDistanceEqualized = real(ifft([HRightDenseSingleDistanceEqualized; ...
    conj(HRightDenseSingleDistanceEqualized(end-1:-1:2, :))]));

%% Select HRIR sparse set (h) for a spherical grid (xl) taken from an icosahedral grid of lower resolution
icoEdgeDivLow = 2;                                      % Subdivision factor of icosahedron's edges
numberOfPointsLow = 10*icoEdgeDivLow^2+2;               % Number of points
fileName = ['icolr', num2str(icoEdgeDivLow), '.mat'];
load(['SphericalGrid\', fileName]);
xl = x*disth(1);
[~, indDir] = max(xh*xl');
errorPositionSelection = max(abs(xh(indDir, :) - xl))
x = xh(indDir, :);
hLeft = hLeftDenseSingleDistanceEqualized(:, indDir);
hRight = hRightDenseSingleDistanceEqualized(:, indDir);

%% Sparse HRTF (H): FFT of sparse HRIR (h)
HLeft = fft(hLeft);     HLeft = HLeft(1:Ns/2+1, :);
HRight = fft(hRight);   HRight = HRight(1:Ns/2+1, :);

%% Spherical Fourier transform of sparse HRTF
L = size(x, 1);                                 % Number of points
N = floor(sqrt(L)-1);                           % Maximum order of spherical harmonics
HLeftnm = sft(HLeft, N, x);
HRightnm = sft(HRight, N, x);

%% Distance varying filter (dn: impulse response, Dn: transfer function)
a = distanceCm*1e-2;                            % Original distance (m)
b = .5;                                         % New distance (m)
regParam = b/a;                                 % Regularization parameter
Dn = zeros(Ns/2+1, (N+1)^2);
for n = 1:N
    for m = -n:n
        Dn(:, n^2+m+n+1) = distanceVaryingFilterSpherical(a, b, n, k, regParam);
    end
end
dn = real(ifft([Dn; conj(Dn(end-1:-1:2, :))]));

%% Distance variation in spherical transform-domain
HLeftNewDistancenm = Dn .* HLeftnm;
HRightNewDistancenm = Dn .* HRightnm;

%% Inverse spherical Fourier transform at original (a) and new (b) distances
azimuthDegree = 90;                         % Desired azimuth in degrees
elevationDegree = 0;                        % Desired elevation in degrees
[aCart(1), aCart(2), aCart(3)] = sph2cart(azimuthDegree*pi/180, elevationDegree*pi/180, a);
[bCart(1), bCart(2), bCart(3)] = sph2cart(azimuthDegree*pi/180, elevationDegree*pi/180, b);
HLefta = isft(HLeftnm, aCart);
HRighta = isft(HRightnm, aCart);
HLeftb = isft(HLeftNewDistancenm, bCart);
HRightb = isft(HRightNewDistancenm, bCart);

%% IFFT of HRTFs (HRIRs)
hLefta = real(ifft([HLefta; conj(HLefta(end-1:-1:2))]));
hRighta = real(ifft([HRighta; conj(HRighta(end-1:-1:2))]));
hLeftb = real(ifft([HLeftb; conj(HLeftb(end-1:-1:2))]));
hRightb = real(ifft([HRightb; conj(HRightb(end-1:-1:2))]));

%% Audio example
[s, FsSignal] = audioread('Audio/SpeechMale48000HzShort.wav');
if FsFilter == FsSignal
    sLefta = fftfilt(hLefta, s);
    sRighta = fftfilt(hRighta, s);
    sLeftb = fftfilt(hLeftb, s);
    sRightb = fftfilt(hRightb, s);
    sBinaurala = [sLefta sRighta];
    sBinauralb = [sLeftb sRightb];
    p = audioplayer([sBinaurala; sBinauralb], FsSignal);
    play(p)
else
    hLeftaResampled = resample(hLefta, FsSignal, FsFilter);
    hRightaResampled = resample(hRighta, FsSignal, FsFilter);
    hLeftbResampled = resample(hLeftb, FsSignal, FsFilter);
    hRightbResampled = resample(hRightb, FsSignal, FsFilter);
    sLefta = fftfilt(hLeftaResampled, s);
    sRighta = fftfilt(hRightaResampled, s);
    sLeftb = fftfilt(hLeftbResampled, s);
    sRightb = fftfilt(hRightbResampled, s);
    sBinaurala = [sLefta sRighta];
    sBinauralb = [sLeftb sRightb];
    p = audioplayer([sBinaurala; sBinauralb], FsSignal);
    play(p)
end

%% Visualization

figure
plot3(xh(:, 1), xh(:, 2), xh(:, 3), 'g.')
hold on
plot3(xl(:, 1), xl(:, 2), xl(:, 3), 'bd')
plot3(x(:, 1), x(:, 2), x(:, 3), 'm*')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
axis equal
axis vis3d
box on
view(viewAngle)

figure
subplot(211)
domain = t*1e3;
plot(domain, equalizerLeft, '.:')
hold on
plot(domain, equalizerRight, '.:')
title('Impulse response $eq(t)$ of diffuse-field equalizer')
xlabel('$t$ (ms)')
ylabel('$eq(t)$')
xlim([domain(1) domain(end)])
legend('Left', 'Right')
grid on
box on

subplot(212)
domain = f*1e-3;
plot(domain, 20*log10(abs(EqualizerLeft)), '.:')
hold on
plot(domain, 20*log10(abs(EqualizerRight)), '.:')
title('Transfer function $Eq(f)$ of diffuse-field equalizer')
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}{|Eq(f)|}$')
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-2:4),'xticklabel', 2.^(-2:4))
xlim([.1 20])
grid on
box on
legend('Left', 'Right')

figure
subplot(211)
domain = t*1e3;
plot(domain, dn, '.:')
title('Impulse response $d_n(t)$ of distance-varying filter')
xlabel('$t$ (ms)')
ylabel('$d_{n}(t)$')
xlim([domain(1) domain(end)])
grid on
box on

subplot(212)
domain = f*1e-3;
value = 20*log10(abs(Dn));
plot(domain, value, '.:')
title('Transfer function $D_n(f)$ of distance-varying filter')
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|D_{n}(f)|$')
xlim([.1 20])
ylim([min(value(:)) max(value(:))])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-1:4), 'xticklabel', 2.^(-1:4))
grid on
box on


figure
subplot(221)
domain = t*1e3;
plot(domain, hLefta)
xlabel('$t$ (ms)')
ylabel('$h(t)$')
xlim([domain(1) domain(end)])
grid on
box on

subplot(222)
domain = t*1e3;
plot(domain, hLeftb)
xlabel('$t$ (ms)')
ylabel('$h(t)$')
xlim([domain(1) domain(end)])
grid on
box on

subplot(223)
domain = f*1e-3;
value = 20*log10(abs(HLefta));
plot(domain, value)
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|H(f)|$')
xlim([.2 20])
ylim([min(value(:)) max(value(:))])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-1:4), 'xticklabel', 2.^(-1:4))
grid on
box on

subplot(224)
domain = f*1e-3;
value = 20*log10(abs(HLeftb));
plot(domain, value)
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|H(f)|$')
xlim([.2 20])
ylim([min(value(:)) max(value(:))])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-1:4), 'xticklabel', 2.^(-1:4))
grid on
box on


% figName = ['figHRIR_azim', num2str(azim(ind)), '_elev', num2str(elev(ind)), '_samples', num2str(Ns)];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

