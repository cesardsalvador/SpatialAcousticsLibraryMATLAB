% 3D Voronoi diagram
% y = voronoi3d(msh, x)
% Input
%  msh : convex hull of x.
%  x   : position vector. size(x) = Nx3.
% Output
%  y   : voronoin diagram. size(y) = 6xM;
% Example
%  msh = voronoin(x);
%  y = voronoi3d(msh, x);
%  plot3(y([1 2], :), y([3 4], :), y([5 6], :),'r-');

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

function y = voronoi3d(msh, x);
tr = TriRep(msh, x(:,1), x(:,2), x(:,3));
numt = size(tr,1);
T = (1:numt)';
neigh = tr.neighbors();
cc = tr.circumcenters();
xcc = cc(:,1);
ycc = cc(:,2);
zcc = cc(:,3);
idx1 = T < neigh(:,1);
idx2 = T < neigh(:,2);
idx3 = T < neigh(:,3);
neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';
y([1 2], :) = xcc(neigh);
y([3 4], :) = ycc(neigh);
y([5 6], :) = zcc(neigh);
