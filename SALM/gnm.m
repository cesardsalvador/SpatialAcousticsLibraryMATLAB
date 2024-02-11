% Spherical spectrum of the Green-Neumann function
% G = gnm(n, m, x, rs, k, kind)
% G = gnm(n, m, x, rs, k, kind, value)
% G = gnm(n, m, x, rs, k, kind, value, norm)
% Inputs:
%  n     : degree associated with the Legendre polynomials (scalar).
%  m     : order in [-n, n] associated with Legendre polynomial (scalar).
%  x     : measurement points x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)]
%  rs    : distance of the sound source
%  k     : wave number
%  kind  : kind of spherical Hankel functions (1: outgoing or 2: incoming)
%  value : 'real' for real-valued (default) or 'complex' for complex-valued.
%  norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% Outputs:
%  G     : spherical expansion coefficients of order n and degree m.
function G = gnm(n, m, x, rs, k, kind, varargin)
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
hnrs = besselhsph(n, kind, k*rs);
dhnr = besselhsph(n-1, kind, k*r(1)) - ...
    (n+1)./(k*r(1)).*besselhsph(n, kind, k*r(1));
f3 = -hnrs./(dhnr.*k*r(1)^2);
G = f1.*f2.*f3;

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
elseif (nargin == 4) & strcmp(varargin{1}, 'none')
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
elseif (nargin == 4) & strcmp(varargin{1}, 'norm')
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
