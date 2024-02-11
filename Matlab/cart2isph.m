% Transform from Cartesian to interaural-polar spherical coordinates.
% The interaural axis is the y-axis. The front position is along the positive x-axis.
%   [pol, lat, dist] = cart2isph(x, y, z)
% Input:
% x, y, z : cartesian coordinates
% Output:
% pol     : polar angle in [-pi, pi] from the xy-plane
% lat     : lateral angle in [-pi/2, pi/2] from the xz-plane
% dist    : radial distance
%
% See also isph2cart.

% Cesar D. Salvador
% cdsalv@gmail.com

function [pol, lat, dist] = cart2isph(x, y, z)
pol = atan2(z, x);
lat = atan2(y, sqrt(x.^2 + z.^2));
dist = sqrt(x.^2 + y.^2 + z.^2);
