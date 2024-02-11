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
