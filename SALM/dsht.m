% Spherical harmonic transform at a new distance (DSHT).
%    Fnm = dsht(F, N, x, d, k, kind)
%    Fnm = dsht(F, N, x, d, k, kind, value)
%    Fnm = dsht(F, N, x, d, k, kind, value, norm)
%    Fnm = dsht(F, N, x, d, k, kind, value, norm, pinv, tol)
%    Fnm = dsht(F, N, x, d, k, kind, value, norm, pinv, tol, quad)
%    [Fnm, G] = dsht(F, N, x, d, k, kind, ...)
%
% Input:
% F     : functions [F1 ... FP] at P positions on the sphere of radius r.
% N     : order of the transform.
% x     : positions x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)]  on the sphere of radius r.
% d     : new distance
% k     : wave number
% kind  : kind of spherical Hankel functions for exp(-jwt) dependence (1: outgoing or 2: incoming)
% value : 'real' for real-valued (default) or 'complex' for complex-valued.
% norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% pinv  : 'moore' for Moore-Penrose pseudo-inverse or 'reg' for regularization.
% tol   : tolerance or regularization parameter to compute the pseudo-inverse
% quad  : quadrature weights for integration on the sphere.
%
% Output:
% Fnm   : spherical spectrum [F00 ... FNM] of order N and length (N+1)^2, at a new distance d.
% G     : matrix of the spherical harmonic transform. Size(G) = (N+1)^2 x P.
%
% See also isht, sht, mpt, impt.

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


function [Fnm, varargout] = dsht(F, N, x, rs, k, kind, varargin)

% Matrix of spherical harmonics at a new distance
H = zeros(size(x, 1), (N+1)^2);
d = zeros(1, (N+1)^2);
if (nargin < 7) || strcmp(varargin{1}, 'real')
    for n = 0:N
        for m = -n:n
            H(:, n^2+n+m+1) = hnm(n, m, x, rs, k, kind, 'real');
            d(n^2+n+m+1) = 1 + n*(n+1);
        end
    end
elseif (nargin == 7) && ( strcmp(varargin{1}, 'real') || ...
                          strcmp(varargin{1}, 'complex') )
    for n = 0:N
        for m = -n:n
            H(:, n^2+n+m+1) = hnm(n, m, x, rs, k, kind, varargin{1});
            d(n^2+n+m+1) = 1 + n*(n+1);
        end
    end
elseif (nargin > 7) && ( strcmp(varargin{2}, 'sch') || ...
                          strcmp(varargin{2}, 'norm') || ...
                          strcmp(varargin{2}, 'none') ) ...
                      && ( strcmp(varargin{1}, 'real') || ...
                          strcmp(varargin{1}, 'complex') )
    for n = 0:N
        for m = -n:n
            H(:, n^2+n+m+1) = hnm(n, m, x, rs, k, kind, varargin{1}, varargin{2});
            d(n^2+n+m+1) = 1 + n*(n+1);
        end
    end
else
    disp('value must be REAL or COMPLEX; and norm must be SCH, NORM or NONE')
    varargout(1) = {[]};
    Fnm = [];
    return
end

% Pseudo-inverse of matrix of spherical harmonics at a new distance
if (nargin < 9)
    varargout(1) = {pinv(H)};
elseif ( (nargin > 8) && (nargin < 11) && (varargin{4} >= 0) && strcmp(varargin{3}, 'moore') )
    varargout(1) = {pinv(H, varargin{4})};
elseif ( (nargin > 8) && (nargin < 11) && (varargin{4} >= 0) && strcmp(varargin{3}, 'reg') )
    varargout(1) = { ( H'*H + varargin{4}*norm(H)*diag(d) ) \ H' };
elseif ( (nargin == 11) && ((sum(varargin{5})-1) < 1e-6) && (varargin{4} >= 0) && strcmp(varargin{3}, 'reg') )
    varargout(1) = { ( H'*diag(varargin{5})*H + varargin{4}*norm(H)*diag(d) ) \ (H'*diag(varargin{5})) };
else
    disp('tol must be positive; pinv must be MOORE or REG; and sum(quad) = 1.')
    varargout(1) = {[]};
    Fnm = [];
end

% Spherical harmonic transform at a new distance
Fnm = (varargout{1}*F.').';



% Functions

% Spherical harmonics at a new distance
% H = Hnm(n, m, x, rs, k, kind)
% H = Hnm(n, m, x, rs, k, kind, value)
% H = Hnm(n, m, x, rs, k, kind, value, norm)
% Inputs:
%  n     : degree associated with the Legendre polynomials (scalar).
%  m     : order in [-n, n] associated with Legendre polynomial (scalar).
%  x     : positions x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)]
%  rs    : new distance
%  k     : wave number
%  kind  : kind of spherical Hankel functions (1: outgoing or 2: incoming)
%  value : 'real' for real-valued (default) or 'complex' for complex-valued.
%  norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% Outputs:
%  H     : spherical harmonics of kind 'kind', degree n, and order m.
function H = hnm(n, m, x, rs, k, kind, varargin)
[azim, elev, r] = cart2sph(x(:,1), x(:,2), x(:,3));
mm = abs(m);
if (nargin < 7) || strcmp(varargin{1}, 'real')
    if m < 0
        f1 = sin(m*azim);
    else
        f1 = cos(m*azim);
    end
elseif (nargin >= 7) && strcmp(varargin{1}, 'complex')
    f1 = exp(1j*m*azim);
end
if (nargin < 8) || strcmp(varargin{2}, 'sch')
    f2 = pnm(n, mm, sin(elev), 'sch');
elseif nargin == 8
    f2 = pnm(n, mm, sin(elev), varargin{2});
end
hnr = besselhsph(n, kind, k*r);
hnrs = besselhsph(n, kind, k*rs);
f3 = hnr./hnrs;
H = f1.*f2.*f3;

% Spherical Hankel function
%        h = besselhsph(n, k, x)
% Input:
%	n     : order
%   kind  : kind (1: outgoing or 2: incoming)
%	x     : argument
% Output
%	h     : spherical Hankel function
function h = besselhsph(m, kind, x)
h = sqrt(pi./(2*x)).*besselh(m+0.5, kind, x);

% Lengendre associated function
% p = pnm(n, m, x)
% p = pnm(n, m, x, nrm)
% Inputs:
%    n : degree (scalar)
%    m : order 0,1,...,n (scalar)
%    x : domain values (array)
% norm : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% Outputs:
%    p : Lengendre associated function of degree n and order m
function p = pnm(n, m, x, varargin)
nn = abs(n);
mm = abs(m);
if (nargin < 4) || strcmp(varargin{1}, 'sch')
    if n == 0
        p = legendre(n, x(:), 'sch');
    elseif n > 0
        Pp = legendre(n, x(:), 'sch');
        if m >= 0
            p = Pp(m+1,:);
        else
            p = (-1)^mm*factorial(n-mm)/factorial(n+mm)*Pp(mm+1,:);
        end
    elseif n < 0
        Pn = legendre(nn-1, x(:), 'sch');
        if m >= 0
            p = Pn(m+1,:);
        else
            p = (-1)^mm*factorial(nn-mm)/factorial(nn+mm)*Pn(mm+1,:);
        end
    end
elseif (nargin == 4) && strcmp(varargin{1}, 'none')
    if n == 0
        p = legendre(n, x(:));
    elseif n > 0
        Pp = legendre(n, x(:));
        if m >= 0
            p = Pp(m+1,:);
        else
            p = (-1)^mm*factorial(n-mm)/factorial(n+mm)*Pp(mm+1,:);
        end
    elseif n < 0
        Pn = legendre(nn-1, x(:));
        if m >= 0
            p = Pn(m+1,:);
        else
            p = (-1)^mm*factorial(nn-mm)/factorial(nn+mm)*Pn(mm+1,:);
        end
    end
elseif (nargin == 4) && strcmp(varargin{1}, 'norm')
    if n == 0
        p = legendre(n, x(:), varargin{1});
    elseif n > 0
        Pp = legendre(n, x(:), varargin{1});
        if m >= 0
            p = Pp(m+1,:);
        else
            p = (-1)^mm*factorial(n-mm)/factorial(n+mm)*Pp(mm+1,:);
        end
    elseif n < 0
        Pn = legendre(nn-1, x(:), varargin{1});
        if m >= 0
            p = Pn(m+1,:);
        else
            p = (-1)^mm*factorial(nn-mm)/factorial(nn+mm)*Pn(mm+1,:);
        end
    end
end
p = reshape(p, size(x));
