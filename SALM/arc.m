% 3D arc between two points passing through a third point
%     T = arc(p1, p2, p0, n)
%
%    INPUT
% p1:   first point [x1, y1, z1]
% p2:   second point [x2, y2, z2]
% p0:   intermediate point [x0, y0, z0]
% n :   number of points on the arc
%    OUTPUT
% T:    arc

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024

function T = arc(p1, p2, p0, n)

% The engine
t = p2-p0; u = p1-p0; v = p1-p2;
w = cross(t,u);
t2 = sum(t.^2); u2 = sum(u.^2); w2 = sum(w.^2);
c = p0+(t2*sum(u.*v)*u-u2*sum(t.*v)*t)/(2*w2); % <-- The center
r = 1/2*sqrt(t2*u2*sum(v.^2)/w2); % <-- The radius
a = p1-c; a = a/norm(a);
b = cross(w,a); b = b/norm(b);
ang = linspace(0,mod(atan2(dot(p2-c,b),dot(p2-c,a)),2*pi),n).';
T = bsxfun(@plus,r*(cos(ang)*a+sin(ang)*b),c);
