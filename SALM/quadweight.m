% Quadrature weights for integration with spherical grids
% q = quadweight(x)
% q = quadweight(x, msh)
%
% Input:
%   x	: spherical grid
%   msh : triangular mesh ( e.g. msh = convhulln(x) )
% Output
%   q   : quadrature weights
%
% See also: ynm, sht, isht.

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


function q = quadweight(x, varargin)
if nargin < 2
    msh = convhulln(x);
else
    msh = varargin{1};
end
tr = TriRep(msh, x);
cc = tr.circumcenters();
cc = max(abs(x(:)))/max(abs(cc(:)))*cc;
[q, centroid] = voronoi_data(length(x), x', length(msh), msh', cc');
end

function [ facearea, centroid ] = voronoi_data ( d_num, d_xyz, face_num, ...
  face_d, v_xyz )

%*****************************************************************************80
%
% VORONOI_DATA computes Voronoi areas and centroids directly.
%
%  Discussion:
%
%    This algorithm works directly from the Delaunay triangulation and the
%    Voronoi vertex locations, without having to determine the explicit
%    description of each Voronoi polygon as a list of Voronoi vertices.
%    This means this routine should be significantly faster than
%    VORONOI_AREAS, at least for larger problems.
%
%
%    The Delaunay triangle has vertices D1, D2, D3, and a Voronoi vertex
%    V which we imagine is inside the triangle, although in fact it need
%    not be.
%
%    Include now the Delaunay vertex averages, D12, D23, and D31, and
%    the Delaunay triangle can be divided into 3 pieces, contributing to
%    the Voronoi polygons associated with D1, D2 and D3 as follows:
%
%      D1 gets D1:V:D31 + D1:D12:V
%      D2 gets D2:V:D12 + D2:D23:V
%      D3 gets D3:V:D23 + D3:D31:V
%
%    The advantage of computing things this way is that it can be done
%    directly from the Delaunay triangle information.  If we must first
%    determine the explicit list of Voronoi vertices that form a particular
%    Voronoi polygon, then compute the area, it is cumbersome and expensive
%    to determine this connectivity information.
%
%    Presumably, if things have been done correctly, the cases where
%    the Voronoi vertex is outside of the Delaunay triangle will correspond
%    to situations where some of the areas come out negative, and 
%    sum to the correct final result.
%
%    The centroids are handled similarly.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    11 May 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer D_NUM, the number of data points and Voronoi polygons.
%
%    Input, real D_XYZ(3,D_NUM), the coodinates of the data points.
%
%    Input, integer FACE_NUM, the number of faces, and Voronoi vertices.
%
%    Input, integer FACE_D(3,FACE_NUM), the data points that form
%    each triangle.
%
%    Input, real V_XYZ(3,FACE_NUM), the coordinates of the Voronoi vertices.
%
%    Output, real AREA(D_NUM), the area of each Voronoi polygon.
%
%    Output, real CENTROID(3,D_NUM), the centroids of each Voronoi polygon.
%
  facearea = zeros ( d_num, 1 );
  centroid = zeros ( 3, d_num );

  r = 1;
%
%  Consider each triangle.
%
  for face = 1 : face_num
%
%  Define some things.
%
    p1 = face_d(1,face);
    p2 = face_d(2,face);
    p3 = face_d(3,face);

    d1 = d_xyz ( 1:3, p1 );
    d2 = d_xyz ( 1:3, p2 );
    d3 = d_xyz ( 1:3, p3 );

    d12 = 0.5 * ( d1 + d2 );
    d23 = 0.5 * ( d2 + d3 );
    d31 = 0.5 * ( d3 + d1 );
    v = v_xyz(1:3,face);
%
%  Force all points to lie on the sphere.
%
    d12 = d12 / norm ( d12 );
    d23 = d23 / norm ( d23 );
    d31 = d31 / norm ( d31 );
    v = v / norm ( v );
%
%  Compute orientation and area of the 6 subtriangles.
%
    o1a = stri_vertices_to_orientation ( d1, v, d31 );
    a1a = stri_vertices_to_area ( r, d1, v, d31 );
    c1a = stri_vertices_to_centroid ( r, d1, v, d31 );

    o1b = stri_vertices_to_orientation ( d1, d12, v );
    a1b = stri_vertices_to_area ( r, d1, d12, v );
    c1b = stri_vertices_to_centroid ( r, d1, d12, v );

    o2a = stri_vertices_to_orientation ( d2, v, d12 );
    a2a = stri_vertices_to_area ( r, d2, v, d12 );
    c2a = stri_vertices_to_centroid ( r, d2, v, d12 );

    o2b = stri_vertices_to_orientation ( d2, d23, v );
    a2b = stri_vertices_to_area ( r, d2, d23, v );
    c2b = stri_vertices_to_centroid ( r, d2, d23, v );

    o3a = stri_vertices_to_orientation ( d3, v, d23 );
    a3a = stri_vertices_to_area ( r, d3, v, d23 );
    c3a = stri_vertices_to_centroid ( r, d3, v, d23 );

    o3b = stri_vertices_to_orientation ( d3, d31, v );
    a3b = stri_vertices_to_area ( r, d3, d31, v );
    c3b = stri_vertices_to_centroid ( r, d3, d31, v );
%
%  Contribute to the Voronoi areas.
%
    facearea(p1) = facearea(p1) + o1a * a1a + o1b * a1b;
    facearea(p2) = facearea(p2) + o2a * a2a + o2b * a2b;
    facearea(p3) = facearea(p3) + o3a * a3a + o3b * a3b;
%
%  Contribute to the Voronoi centroids.
%
    centroid(1:3,p1) = centroid(1:3,p1) + o1a * a1a * c1a(1:3) + o1b * a1b * c1b(1:3);
    centroid(1:3,p2) = centroid(1:3,p2) + o2a * a2a * c2a(1:3) + o2b * a2b * c2b(1:3);
    centroid(1:3,p3) = centroid(1:3,p3) + o3a * a3a * c3a(1:3) + o3b * a3b * c3b(1:3);

  end
  
%
%  Scale the areas of the faces
%
   facearea = abs(facearea)*4*pi*max(abs(d_xyz(:)))^2/sum(abs(facearea));
  
%
%  Normalize the centroids to lie on the sphere.
%
  for j = 1 : d_num
    centroid(1:3,j) = centroid(1:3,j) / norm ( centroid(1:3,j) );
  end

  return
end
%
%*****************************************************************************80

function o = stri_vertices_to_orientation ( a, b, c )

%*****************************************************************************80
%
%% STRI_VERTICES_TO_ORIENTATION seeks the orientation of a spherical triangle.
%
%  Discussion:
%
%    Three points on a sphere actually compute two triangles; typically
%    we are interested in the smaller of the two.
%
%    As long as our triangle is "small", we can define an orientation
%    by comparing the direction of the centroid against the normal
%    vector (C-B) x (A-B).  If the dot product of these vectors
%    is positive, we say the triangle has positive orientation.
%
%    By using information from the triangle orientation, we can correctly
%    determine the area of a Voronoi polygon by summing up the pieces
%    of Delaunay triangles, even in the case when the Voronoi vertex
%    lies outside the Delaunay triangle.  In that case, the areas of
%    some of the Delaunay triangle pieces must be formally negative.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    11 May 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A(3), B(3), C(3), three points on a sphere.
%
%    Output, integer O, is +1 if the spherical triangle is judged to
%    have positive orientation, and -1 otherwise.
%

%
%  Please, column vectors only!
%
  a = a( : );
  b = b( : );
  c = c( : );
%
%  Centroid.
%
  cd = ( a + b + c ) / 3.0;
%
%  Cross product ( C - B ) x ( A - B );
%
  v1 = c - b;
  v2 = a - b;

  cp = zeros ( 3, 1 );
  cp(1) = v1(2) * v2(3) - v1(3) * v2(2);
  cp(2) = v1(3) * v2(1) - v1(1) * v2(3);
  cp(3) = v1(1) * v2(2) - v1(2) * v2(1);
%
%  Compare the directions.
%
  if ( cp' * cd < 0 )
    o = - 1;
  else
    o = + 1;
  end

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function area = stri_vertices_to_area ( r, v1, v2, v3 )

%*****************************************************************************80
%
%% STRI_VERTICES_TO_AREA computes the area of a spherical triangle.
%
%  Discussion:
%
%    A sphere centered at 0 in 3D satisfies the equation:
%
%      X * X + Y * Y + Z * Z = R * R
%
%    A spherical triangle is specified by three points on the surface
%    of the sphere.
%
%    The area formula is known as Girard's formula.
%
%    The area of a spherical triangle is:
%
%      AREA = ( A + B + C - PI ) * R * R
%
%    where A, B and C are the (surface) angles of the triangle.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 April 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real R, the radius of the sphere.
%
%    Input, real V1(3), V2(3), V3(3), the vertices of the triangle.
%
%    Output, real AREA, the area of the spherical triangle.
%

%
%  Compute the lengths of the sides of the spherical triangle.
%
  [ as, bs, cs ] = stri_vertices_to_sides ( r, v1, v2, v3 );
%
%  Get the spherical angles.
%
  [ a, b, c ] = stri_sides_to_angles ( r, as, bs, cs );
%
%  Get the area.
%
  area = stri_angles_to_area ( r, a, b, c );

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ as, bs, cs ] = stri_vertices_to_sides ( r, v1, v2, v3 )

%*****************************************************************************80
%
%% STRI_VERTICES_TO_SIDES computes spherical triangle sides.
%
%  Discussion:
%
%    We can use the ACOS system call here, but the ARC_COSINE routine
%    will automatically take care of cases where the input argument is
%    (usually slightly) out of bounds.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 April 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real R, the radius of the sphere.
%
%    Input, real V1(3), V2(3), V3(3), the vertices of the spherical
%    triangle.
%
%    Output, real AS, BS, CS, the (geodesic) length of the sides
%    of the triangle.
%
  dim_num = 3;

  as = r * arc_cosine ( v2(1:dim_num)' * v3(1:dim_num) / r.^2 );
  bs = r * arc_cosine ( v3(1:dim_num)' * v1(1:dim_num) / r.^2 );
  cs = r * arc_cosine ( v1(1:dim_num)' * v2(1:dim_num) / r.^2 );

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = arc_cosine ( c )

%*******************************************************************************
%
%% ARC_COSINE computes the arc cosine function, with argument truncation.
%
%  Discussion:
%
%    If you call your system ACOS routine with an input argument that is
%    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
%    surprise (I did).
%
%    This routine simply truncates arguments outside the range.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 January 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real C, the argument.
%
%    Output, real VALUE, an angle whose cosine is C.
%
  c2 = c;
  c2 = max ( c2, -1.0 );
  c2 = min ( c2, +1.0 );

  value = acos ( c2 );

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ a, b, c ] = stri_sides_to_angles ( r, as, bs, cs )

%*****************************************************************************80
%
%% STRI_SIDES_TO_ANGLES computes spherical triangle angles.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real R, the radius of the sphere.
%
%    Input, real AS, BS, CS, the (geodesic) length of the 
%    sides of the triangle.
%
%    Output, real A, B, C, the spherical angles of the triangle.
%    Angle A is opposite the side of length AS, and so on.
%
  asu = as / r;
  bsu = bs / r;
  csu = cs / r;
  ssu = ( asu + bsu + csu ) / 2.0;

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - asu )     ) );

  a = 2.0 * atan ( tan_a2 );

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) );

  b = 2.0 * atan ( tan_b2 );

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - csu )     ) );

  c = 2.0 * atan ( tan_c2 );

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function area = stri_angles_to_area ( r, a, b, c )

%*****************************************************************************80
%
%% STRI_ANGLES_TO_AREA computes the area of a spherical triangle.
%
%  Discussion:
%
%    A sphere centered at 0 in 3D satisfies the equation:
%
%      X*X + Y*Y + Z*Z = R*R
%
%    A spherical triangle is specified by three points on the surface
%    of the sphere.
%
%    The area formula is known as Girard's formula.
%
%    The area of a spherical triangle is:
%
%      AREA = ( A + B + C - PI ) * R*R
%
%    where A, B and C are the (surface) angles of the triangle.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real R, the radius of the sphere.
%
%    Input, real A, B, C, the angles of the triangle.
%
%    Output, real AREA, the area of the spherical triangle.
%

%
%  Apply Girard's formula.
%
  area = r * r * ( a + b + c - pi );

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vs = stri_vertices_to_centroid ( r, v1, v2, v3 )

%*****************************************************************************80
%
%% STRI_VERTICES_TO_CENTROID_3D gets a spherical triangle centroid in 3D.
%
%  Discussion:
%
%    A sphere centered at 0 in 3D satisfies the equation:
%
%      X*X + Y*Y + Z*Z = R*R
%
%    A spherical triangle is specified by three points on the sphere.
%
%    The (true) centroid of a spherical triangle is the point
%
%      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
%
%    Note that the true centroid does NOT, in general, lie on the sphere.  
%
%    The "flat" centroid VF is the centroid of the planar triangle defined by
%    the vertices of the spherical triangle.
%
%    The "spherical" centroid VS of a spherical triangle is computed by
%    the intersection of the geodesic bisectors of the triangle angles.
%    The spherical centroid lies on the sphere.
%
%    VF, VT and VS lie on a line through the center of the sphere.  We can
%    easily calculate VF by averaging the vertices, and from this determine
%    VS by normalizing.
%
%    Of course, we still will not have actually computed VT, which lies
%    somewhere between VF and VS!
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real R, the radius of the sphere.
%
%    Input, real V1(3), V2(3), V3(3), the vertices of the triangle.
%
%    Output, real VS(3), the coordinates of the "spherical
%    centroid" of the spherical triangle.
%
  dim_num = 3;

  vs(1:dim_num,1) = ( v1(1:dim_num) + v2(1:dim_num) + v3(1:dim_num) ) / 3.0;

  norm = sqrt ( sum ( vs(1:dim_num,1).^2 ) );

  vs(1:dim_num,1) = r * vs(1:dim_num,1) / norm;

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

