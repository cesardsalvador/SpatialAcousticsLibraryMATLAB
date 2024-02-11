% Transform spherical coordinates to ISO spherical coordinates
%   [r, theta, phi] = sph2sphiso(THETA, PHI, R)
% Input:
% THETA   : azimuth in [0, 2*pi]
% PHI     : elevation in [-pi/2, pi/2]
% R       : radius
% Output:
% r       : radius
% theta   : inclinaton or zenith in [0, pi]
% phi     : azimuth in [0, 2*pi]

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [r, theta, phi] = sph2sphiso(THETA, PHI, R)
phi = THETA; ind = find(phi<0); phi(ind) = 2*pi+phi(ind);
theta = pi/2 - PHI;
r = R;
