% Transform Cartesian to interaural spherical coordinates.
% The interaural axis is the y-axis. The front position is along the positive x-axis.
%   [alpha, beta, r] = cart2sphiso_old(x, y, z)
% Input:
% x, y, z : cartesian coordinates
% Output:
% r       : radius
% alpha   : azimuthal or lateral angle in [-pi/2, pi/2]
% beta    : elevation or polar angle in [-pi, pi]

function [alpha, beta, r] = cart2isph_old(x, y, z)
alpha = atan2(y, sqrt(x.^2 + z.^2));
beta = atan2(z, x);
r = sqrt(x.^2 + y.^2 + z.^2);
