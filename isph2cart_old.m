% Transform interaural spherical coordinates to Cartesian coordinates.
% The interaural axis is the y-axis. The front position is along the positive x-axis.
%   [x, y, z] = isph2cart_old(alpha, beta, r)
% Input:
% alpha   : azimuthal or lateral angle in [-pi/2, pi/2]
% beta    : elevation or polar angle in [-pi, pi]
% r       : radial distance
% Output:
% x, y, z : cartesian coordinates

function [x, y, z] = isph2cart_old(alpha, beta, r)
x = r.*cos(alpha).*cos(beta);
y = r.*sin(alpha);
z = r.*cos(alpha).*sin(beta);
