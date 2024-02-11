% Transform Cartesian to spherical coordinates (ISO)
%   [r, theta, phi] = cart2sphiso(x, y, z)
% Input:
% x, y, z : cartesian coordinates
% Output:
% r       : radius
% theta   : inclinaton or zenith in [0, pi]
% phi     : azimuth in [0, 2*pi]

function [r, theta, phi] = cart2sphiso(x, y, z)
r = sqrt(x.^2 + y.^2 + z.^2);
theta = atan2(sqrt(x.^2 + y.^2), z);
phi = atan2(y, x);
ind = phi<0;
phi(ind) = 2*pi+phi(ind);