% Random distribution of points inside a sphere (randsph).
% Uniform distribution along angles to avoid concentration around the poles.
% Quadratic distribution along distance to avoid concentration around the origin.
%      x = randsph(N, r, type)
%
%   Input:
%   N     : number of points.
%   r     : radius of the sphere
%   type  : 'surface' or 'volume'
%
%   Output:
%   x     : N positions in Cartesian coordinates.

% June 4, 2020
% Cesar D. Salvador
% daniel@si-in.com

% Reference
% Weisstein, Eric W. "Sphere Point Picking."
% From MathWorld--A Wolfram Web Resource
% https://mathworld.wolfram.com/SpherePointPicking.html

function x = randsph(N, r, type)
u = randn(N, 3);
u = diag( 1./sqrt(u(:, 1).^2 + u(:, 2).^2 + u(:, 3).^2) ) * u;
switch type
    case 'surface'
        x = r*u;
    case 'volume'
        r0 = r - ( 2.^(log2(r/N)+(log2(r)-log2(r/N))*rand(N, 1)) );
        x = diag(r0) * u;
end
end
