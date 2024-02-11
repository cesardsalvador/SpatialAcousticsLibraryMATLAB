%% Head-related impulse response (HRIR) from IRCAM Listen database
% This script reads and plots the HRIR downloaded from
% http://recherche.ircam.fr/equipes/salles/listen/download.html

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


%% Initialization
close all
clearvars
clc

%% Format for vizualisation
cmap = jet(256);
font_size = 18;
set(0, ...
    'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmap, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', font_size, ...
    'defaultAxesFontSize', font_size, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')
viewAngle = [-40 10];

%% Format for saving figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Load HRIR (hLefta and hRighta)
subjectNumber = 53;                 % subject number (reduced set: 8, 13, 22, 31, 32, 48, 53)
typeHRIR = 'COMPENSATED';           % 'COMPENSATED' or 'RAW'
headModelName = ['listenIRCAMsubject', sprintf('%02d', subjectNumber)];
disth = 1.95;                       % radial distance in m
xEarLeft = [0 0.085 0];             % Left-ear position
xEarRight = [0 -0.085 0];           % Right-ear position
sampleCircshift = [0 0];

folderName = ['C:\Prog\Matlab\Data\hrir\ircam\SUBJECTS\IRC_10', sprintf('%02d', subjectNumber), '\', typeHRIR, '\MAT\HRIR\'];
fileName = ['IRC_10', sprintf('%02d', subjectNumber), '_', typeHRIR(1),'_HRIR.mat'];
load([folderName, fileName])
switch typeHRIR
    case 'RAW'
        [xh(:, 1), xh(:, 2), xh(:, 3)] = sph2cart( ...      % Positions of sound sources
            l_hrir_S.azim_v*pi/180, ...
            l_hrir_S.elev_v*pi/180, disth);
        hLeft = l_hrir_S.content_m';                       % Left-ear HRIR
        hRight = r_hrir_S.content_m';                      % Right-ear HRIR
        Ns = size(hLeft, 1);                               % Number of samples
        FsFilter = l_hrir_S.sampling_hz;                         % Sampling frequency
    case 'COMPENSATED'
        [xh(:, 1), xh(:, 2), xh(:, 3)] = sph2cart( ...     % Positions of sound sources
            l_eq_hrir_S.azim_v*pi/180, ...
            l_eq_hrir_S.elev_v*pi/180, disth);
        hLeft = l_eq_hrir_S.content_m';                    % Left-ear HRIR
        hRight = r_eq_hrir_S.content_m';                   % Right-ear HRIR
        Ns = size(hLeft, 1);                               % Number of samples
        FsFilter = l_eq_hrir_S.sampling_hz;                      % Sampling frequency
end

%% One-dimensional domains
t = (0:Ns-1)' * 1/FsFilter;                                       % Time domain
f = (0:Ns/2)' * FsFilter/Ns;                                      % Frequency domain
c = 344;                                                    % Speed of sound in air (m/s)
k = 2*pi*f/c;                                               % Wave number

%% Select HRIRs for two directions
azim1 = 90;
elev1 = 0;
azim2 = -90;
elev2 = 0;
[x1(1), x1(2), x1(3)] = sph2cart(azim1*pi/180, elev1*pi/180, disth);
[x2(1), x2(2), x2(3)] = sph2cart(azim2*pi/180, elev2*pi/180, disth);
[~, indDir1] = max(xh*x1');
[~, indDir2] = max(xh*x2');
hLeft1 = hLeft(:, indDir1);
hRight1 = hRight(:, indDir1);
hLeft2 = hLeft(:, indDir2);
hRight2 = hRight(:, indDir2);

%% Read monoaural audio files
[s1, FsSignal] = audioread('Audio/SpeechMale48000Hz.wav');
[s2, ~] = audioread('Audio/SpeechFemale48000Hz.wav');
hLeft1 = resample(hLeft1, FsSignal, FsFilter);
hRight1 = resample(hRight1, FsSignal, FsFilter);
hLeft2 = resample(hLeft2, FsSignal, FsFilter);
hRight2 = resample(hRight2, FsSignal, FsFilter);
s1Left = fftfilt(hLeft1, s1);
s1Right = fftfilt(hRight1, s1);
s2Left = fftfilt(hLeft2, s2);
s2Right = fftfilt(hRight2, s2);
sBinaural = [s1Left+s2Left s2Left+s2Right];
sBinaural = sBinaural/(sqrt(2)*max(abs(sBinaural(:))));

%% Play binaural audio
p = audioplayer(sBinaural, FsSignal);
play(p)

%% Write binaural audio file
%audiowrite('Audio/binauralDuoSpeech48000Hz.wav', sBinaural, FsSignal);

%% Visualization
plot3(xh(:, 1), xh(:, 2), xh(:, 3), 'g.')
hold on
plot3(x1(1), x1(2), x1(3), 'm*')
plot3(x2(1), x2(2), x2(3), 'bo')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
axis equal
axis vis3d
box on
view(viewAngle)

figure
domain = t*1e3;
plot(domain, hLeft, 'r')
hold on
plot(domain, hRight, 'b')
title(['Fs=', num2str(FsFilter), ' Hz, Samples= ', num2str(Ns)])
xlabel('$t$ (ms)')
ylabel('Amplitude')
xlim([domain(1) domain(end)])
maxLim = max(abs([hLeft(:); hRight(:)]));
ylim([-maxLim maxLim])


% fig_name = ['figHRIR_azim', num2str(azim(ind)), '_elev', num2str(elev(ind)), '_samples', num2str(Ns)];
% print(figFormat, figRes, [figFolder, fig_name, figNameExt])