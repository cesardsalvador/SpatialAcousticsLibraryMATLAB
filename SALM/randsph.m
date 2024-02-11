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

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

% Reference and citation
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.
%     DOI: 10.1109/TASLP.2017.2778562
% [2] C. D. Salvador et al., “Design theory for binaural synthesis:
%     Combining microphone array recordings and head-related transfer function datasets,”
%     Acoust. Sci. Technol., vol. 38, no. 2, pp. 51–62, Mar. 2017.
%     DOI: 10.1250/ast.38.51


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
