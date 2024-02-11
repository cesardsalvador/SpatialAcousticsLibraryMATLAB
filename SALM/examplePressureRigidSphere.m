
close all
clearvars
clc

% Parameters for vizualisation
%cmap = othercolor('BuDRd_18', 256);
cmap = jet(256);
fontSize = 18;
set(0, 'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmap, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', fontSize, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')

fMin = 100;
fMax = 20e3;
fTick = 5:5:15;
fTick_log = 2.^(-1:4);

azimMic = 30;
elevMic = 30;
radiusMic = 0.085;
[xMic(:, 1), xMic(:, 2), xMic(:, 3)] = sph2cart(azimMic*pi/180, elevMic*pi/180, radiusMic);

azimSource = 0;
elevSource = 0;
radiusSource = 1.5;
[xSource(:, 1), xSource(:, 2), xSource(:, 3)] = sph2cart(azimSource*pi/180, elevSource*pi/180, radiusSource);

% Angles between sound source positions and microphones positions
cosTheta = xMic*xSource'./(sqrt(sum(xMic.^2, 2))*sqrt(sum(xSource.^2, 2))');
cosTheta(abs(cosTheta)>=1) = sign(cosTheta(abs(cosTheta)>=1));
Theta = acosd(cosTheta);                                           % 3D angle

Ns = 512;
Fs = 48e3;
c = 343;
f = (0:Ns/2)*Fs/Ns;
k = 2*pi*f/c;
t = (0:Ns-1)*1/Fs;
thres = 1e-4;

PSphericalWave = zeros(Ns/2+1, 1);
PPlaneWave = zeros(Ns/2+1, 1);
for I = 2:Ns/2+1
    PSphericalWave(I) = conj(pressureRigidSphere(radiusMic, radiusSource, Theta, k(I), thres));
%    PPlaneWave(I) = conj(pressureRigidSpherePlaneWave(radiusMic, Theta, k(I), thres));
    PPlaneWave(I) = pressureRigidSpherePlaneWave(radiusMic, Theta, k(I), thres);
end
pSphericalWave = real(ifft([ PSphericalWave; conj(PSphericalWave(end-1:-1:2)) ]));
pPlaneWave = real(ifft([ PPlaneWave; conj(PPlaneWave(end-1:-1:2)) ]));

figure()
subplot(211)
plot(f*1e-3, 20*log10(abs(PSphericalWave)))
xlabel('$f$ (kHz)')
ylabel('$20 \log_{10}|P(f)|$')
xlim([fMin fMax]*1e-3)
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', fTick_log, 'xticklabel', fTick_log)
box on

subplot(212)
plot(t*1e3, pSphericalWave)
xlabel('$t$ (ms)')
ylabel('$p(t)$')
set(gca, 'xscale', 'linear', 'yscale', 'linear')
axis tight
box on


figure
subplot(211)
plot(f*1e-3, 20*log10(abs(PPlaneWave)))
xlabel('$f$ (kHz)')
ylabel('$20 \log_{10}|P(f)|$')
xlim([fMin fMax]*1e-3)
set(gca, 'xscale', 'log', 'yscale', 'linear')
set(gca, 'xtick', fTick_log, 'xticklabel', fTick_log)
box on

subplot(212)
plot(t*1e3, pPlaneWave)
xlabel('$t$ (ms)')
ylabel('$p(t)$')
set(gca, 'xscale', 'linear', 'yscale', 'linear')
axis tight
box on
