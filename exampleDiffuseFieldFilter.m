%% Example of diffuse-field equalization
%
% Summary:
%   1) Reads a dense spherical dataset of HRIRs from a SOFA file.
%   2) Applies diffuse field equalization to the HRIR dense dataset.
%
% Notation:
%   h: head-related impulse response (HRIR)
%   H: head-related transfer function (HRTF)
%
% Requirements:
%   1) SOFA Matlab/Octave API version 1.0.2
%      https://sourceforge.net/projects/sofacoustics/
%   2) Near-distance HRTF dataset
%      https://cesardsalvador.github.io/download.html%
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
% [2] J. C. Middlebrooks, “Individual differences in external-ear transfer functions reduced
% by scaling in frequency,” J. Acoust. Soc. Am., vol. 106, no. 3, pp. 1480–1492, 1999.
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
sourceDistribution = 'sph_ico_dir2562_dist21';
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

%% Cirular shift of HRIR dense set
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

%% Select HRIR dense spherical grid (xh) at a single distance
distanceCm = 200;                                 % Distance (cm) from listener to point source
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
weightDense = quadweight(xhUnitSphere)'/(4*pi);

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

%% Audio example
azimDegree = 90;                        % Desired azimuth in degrees
elevDegree = 10;                        % Desired elevation in degrees
[y(:, 1), y(:, 2), y(:, 3)] = sph2cart(azimDegree*pi/180, elevDegree*pi/180, distanceCm*1e-2);
[~, indPosition] = max(abs(xh*y'));
[s, FsSignal] = audioread('Audio/SpeechMale48000HzShort.wav');
if FsFilter == FsSignal
    sLeft = fftfilt(hLeftDenseSingleDistance(:, indPosition), s);
    sRight = fftfilt(hRightDenseSingleDistance(:, indPosition), s);
    sLeftEqualized = fftfilt(hLeftDenseSingleDistanceEqualized(:, indPosition), s);
    sRightEqualized = fftfilt(hRightDenseSingleDistanceEqualized(:, indPosition), s);
    sBinaural = [sLeft sRight];
    sBinauralEqualized = [sLeftEqualized sRightEqualized];
    p = audioplayer([sBinaural; sBinauralEqualized], FsSignal);
    play(p)
else
    hLeftResampled = resample(hLeftDenseSingleDistance(:, indPosition), FsSignal, FsFilter);
    hRightResampled = resample(hRightDenseSingleDistance(:, indPosition), FsSignal, FsFilter);
    hLeftResampledEqualized = resample(hLeftDenseSingleDistanceEqualized(:, indPosition), FsSignal, FsFilter);
    hRightResampledEqualized = resample(hRightDenseSingleDistanceEqualized(:, indPosition), FsSignal, FsFilter);
    sLeft = fftfilt(hLeftResampled, s);
    sRight = fftfilt(hRightResampled, s);
    sLeftEqualized = fftfilt(hLeftResampledEqualized(:, indPosition), s);
    sRightEqualized = fftfilt(hRightResampledEqualized(:, indPosition), s);
    sBinaural = [sLeft sRight];
    sBinauralEqualized = [sLeftEqualized sRightEqualized];
    p = audioplayer([sBinaural; sBinauralEqualized], FsSignal);
    play(p)
end

%% Visualization

figure
plot3(xh(:, 1), xh(:, 2), xh(:, 3), 'g.')
hold on
plot3(y(1), y(2), y(3), 'rd')
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
legend('Left', 'Right', 'Location', 'southwest')


figure
subplot(221)
domain = t*1e3;
plot(domain, hLeftDenseSingleDistance)
xlabel('$t$ (ms)')
ylabel('$h(t)$')
xlim([domain(1) domain(end)])

subplot(222)
domain = t*1e3;
plot(domain, hLeftDenseSingleDistanceEqualized)
xlabel('$t$ (ms)')
ylabel('$h_{\rm{Equalized}}(t)$')
xlim([domain(1) domain(end)])

subplot(223)
domain = f*1e-3;
value = 20*log10(abs(HLeftDenseSingleDistance));
plot(domain, value)
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|H(f)|$')
xlim([.1 20])
ylim([min(value(:)) max(value(:))])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-1:4), 'xticklabel', 2.^(-1:4))
set(gca, 'ytick', -120:24:120, 'yticklabel', -120:24:120)
grid on
box on

subplot(224)
domain = f*1e-3;
value = 20*log10(abs(HLeftDenseSingleDistanceEqualized));
plot(domain, value)
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|H_{\rm{Equalized}}(f)|$')
xlim([.1 20])
ylim([min(value(:)) max(value(:))])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-1:4), 'xticklabel', 2.^(-1:4))
set(gca, 'ytick', -120:24:120, 'yticklabel', -120:24:120)
grid on
box on

% figName = ['figHRIR_azim', num2str(azim(ind)), '_elev', num2str(elev(ind)), '_samples', num2str(Ns)];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

