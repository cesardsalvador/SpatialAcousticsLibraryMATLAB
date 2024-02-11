% Fourier-Legendre transform (FLT).
%    Fn = flt(F, N, theta)
%    Fn = flt(F, N, theta, norm)
%    Fn = flt(F, N, theta, norm, pinv, tol)
%    [Fn, G] = sht(F, N, theta, ...)
%
% Input:
% F     : vector function F = [F1 ... FM] evaluated at M angles.
% N     : order of the transform (natural number).
% theta : angles theta = [theta1 theta2 ... thetaM].
% norm : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% pinv  : 'moore' for Moore-Penrose pseudo-inverse (default), or 'reg' for regularization.
% tol   : tolerance or regularization parameter to compute the pseudo-inverse.
%
% Output:
% Fn    : Fourier-Legendre transform [F_0 ... F_(N-1)] of order N.
% G     : matrix of the Fourier-Legendre transform. Size(G) = N x M.
%
% See also iflt, sht, isht.

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [Fn, varargout] = flt(F, N, theta, varargin)

% Matrix of Legendre polynomials
H = zeros(length(theta), N);
if (nargin < 4)
    for n = 0:N-1
            H(:, n+1) = pnm(n, 0, sin(theta));
    end
elseif (nargin >= 4) && ( strcmp(varargin{1}, 'sch') || ...
                          strcmp(varargin{1}, 'norm') || ...
                          strcmp(varargin{1}, 'none') )
    for n = 0:N-1
            H(:, n+1) = pnm(n, 0, sin(theta), varargin{1});
    end
else
    disp('norm must be SCH, NORM or NONE')
    varargout(1) = {[]};
    Fn = [];
    return
end

% Pseudo-inverse of matrix of Legendre polynomials
if (nargin < 5)
    varargout(1) = {pinv(H)};
elseif ( (nargin > 4) && (nargin < 7) && (varargin{3} >= 0) && strcmp(varargin{2}, 'moore') )
    varargout(1) = {pinv(H, varargin{3})};
elseif ( (nargin > 4) && (nargin < 7) && (varargin{3} >= 0) && strcmp(varargin{2}, 'reg') )
    varargout(1) = { ( H'*H + varargin{3}*norm(H) ) \ H' };
else
    disp('tol must be positive; and pinv must be MOORE or REG.')
    varargout(1) = {[]};
    Fn = [];
end

% Fourier-Legendre transform
Fn = (varargout{1}*F.').';


