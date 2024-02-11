% Transform spherical coordinates to ISO spherical coordinates
%   [THETA, PHI, R] = sphiso2sph(r, theta, phi)
% Input:
% r       : radius
% theta   : inclinaton or zenith in [0, pi]
% phi     : azimuth in [0, 2*pi]
% Output:
% THETA   : azimuth in [0, 2*pi]
% PHI     : elevation in [-pi/2, pi/2]
% R       : radius

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [THETA, PHI, R] = sphiso2sph(r, theta, phi)
THETA = phi; ind = find(THETA<0); THETA(ind) = pi-THETA(ind);
PHI = pi/2 - theta;
R = r;
