%% Example of feedbackDelayImpulseResponse
% This script designs and plots an artificial room impulse response as the
% impulse response of a feedback delay network

% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/
% December 15, 2020

%% Clean and prepare the Matlab environment
close all                                   % Close all figures
clearvars                                   % Clear all variables in the workspace
clc                                         % Clean the command window
workspace                                   % Make sure the workspace panel is showing.
format longg
format compact

%% Format for vizualisation
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

%% Format for saving figures
figFolder = 'Fig\';
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
%figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';

%% Parameters
T60 = 0.7;                                              % Reverberation time in seconds
DRRdB = 20;                                             % Direct-to-reverberant energy ratio in dB
FsFilter = 48e3;                                        % Sampling frequency
Ns = ceil(T60*FsFilter);                                % Number of samples
L = 8;                                                  % Number of delay lines
minDelayMiliSecond = 10;                                % Minimum delay in miliseconds
maxDelayMiliSecond = 40;                                % Maximum delay in miliseconds

%% One-dimensional domains
t = (0:Ns-1)' * 1/FsFilter;                                   % Time domain
f = (0:Ns/2)' * FsFilter/Ns;                                  % Frequency domain

%% Vectors and matrices to design the feedback delay network
minDelaySample = floor(minDelayMiliSecond*1e-3*FsFilter);     % Minimum delay in samples
maxDelaySample = floor(maxDelayMiliSecond*1e-3*FsFilter);     % Maximum delay in samples
delay = randi([minDelaySample maxDelaySample], [L 1]);
attenuationPerSampledB = -60/(T60*FsFilter);
attenuationPerSample = 10^(attenuationPerSampledB/20);
attenuation = attenuationPerSample.^delay;
attenuationMatrix = diag(attenuation);
scalarFeedbackMatrix = hadamard(L)/sqrt(L);
feedbackMatrix = scalarFeedbackMatrix * attenuationMatrix;
inputGain = ones(L, 1);
directGain = 1;
DRR = 10^(-DRRdB/20);
outputGain = DRR * attenuation/max(attenuation);

%% Room impulse response (RIR) generated with a feedback delay network
RIR = feedbackDelayImpulseResponse(delay, feedbackMatrix, inputGain, directGain, outputGain, Ns, FsFilter);
RTF = fft(RIR); RTF = RTF(1:Ns/2+1);

%% Audio example
[s, FsSignal] = audioread('Audio/SpeechFemale48000HzShort.wav');
RIRresampled = resample(RIR, FsSignal, FsFilter);
sReverb = fftfilt(RIRresampled, s);
sExample = [s; sReverb];
sExample = sExample / max(abs(sExample(:)));
p = audioplayer(sExample, FsSignal);
play(p)

%% Visualization
figure
subplot(311)
domain = t*1e3;                                     % Time domain in ms (x axis)
value = 20*log10(abs(RIR));                          % Impulse response in dB (y axis)
plot(domain, value, '.:')
title('Room impulse response')
xlabel('$t$ (ms)')
ylabel('$20\log_{10}|h(t)|$')
xlim([min(domain) max(domain)])
ylim([-DRRdB-60 max(value)])
set(gca, 'xscale', 'linear', 'yscale', 'linear')
set(gca, 'ytick', -120:20:60, 'yticklabel', -120:20:60)
grid on
box on

subplot(312)
domain = f*1e-3;                                    % Frequency domain in kHz (x axis)
value = 20*log10(abs(RTF));                          % Magnitude of transfer function in dB (y axis)
plot(domain, value)
title('Magnitude of room transfer function')
xlabel('$f$ (kHz)')
ylabel('$20\log_{10}|H(f)|$')
xlim([.1 20])
ylim([min(value) max(value)])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
%set(gca, 'ytick', -24:6:24, 'yticklabel', -24:6:24)
grid on
box on

subplot(313)
domain = f(2:end)*1e-3;                                 % Frequency domain in kHz (x axis)
value = -1e3 * diff(unwrap(angle(RTF)))/(2*pi*FsFilter/Ns);	% Group delay in miliseconds (y axis)
plot(domain, value)
title('Group delay of room transfer function')
xlabel('$f$ (kHz)')
ylabel('$-\frac{\partial{\Phi(\omega)}}{\partial{\omega}}$ (ms)')
xlim([.1 20])
ylim([min(value) max(value)])
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
%set(gca, 'ytick', -120:10:60, 'yticklabel', -120:10:60)
grid on
box on

%% Saving figures
% figName = ['figFourierTransform', num2str(f0), 'HzTone'];
% print(figFormat, figRes, [figFolder, figName, figNameExt])

