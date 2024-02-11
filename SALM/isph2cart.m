% Transform from interaural-polar spherical coordinates to Cartesian coordinates.
% The interaural axis is the y-axis. The front position is along the positive x-axis.
%   [x, y, z] = isph2cart(pol, lat, dist)
% Input:
% pol     : polar angle in [-pi, pi] from the xy-plane
% lat     : lateral angle in [-pi/2, pi/2] from the xz-plane
% dist    : radial distance
% Output:
% x, y, z : cartesian coordinates
%
% See also cart2isph.

% Cesar D. Salvador
% cdsalv@gmail.com

function [x, y, z] = isph2cart(pol, lat, dist)
x = dist.*cos(lat).*cos(pol);
y = dist.*sin(lat);
z = dist.*cos(lat).*sin(pol);
