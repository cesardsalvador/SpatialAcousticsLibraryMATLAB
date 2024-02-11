% Transform ISO spherical coordinates to Cartesian coordinates
%   [x, y, z] = sphiso2cart(r, theta, phi)
% Input:
% r       : radius
% theta   : inclinaton or zenith in [0, pi]
% phi     : azimuth in [0, 2*pi]
% Output:
% x, y, z : cartesian coordinates

function [x, y, z] = sphiso2cart(r, theta, phi)
x = r.*cos(phi).*sin(theta);
y = r.*sin(phi).*sin(theta);
z = r.*cos(theta);
