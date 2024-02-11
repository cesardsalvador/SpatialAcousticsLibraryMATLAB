% Inverse Fourier Legendre transform (IFLT).
%    F = iflt(Fn, theta)
%    F = iflt(Fn, theta, norm)
%    F = iflt(Fn, theta, norm, L)
%    [F, H] = iflt(Fn, theta, ...)
%
% Input:
% Fn    : Fourier-Legendre transform [F_0 ... F_(N-1)] of order N.
% theta : angles theta = [theta1 theta2 ... thetaM].
% norm  : 'sch' (default Schmidt seminorm), 'norm' (full norm), or 'none' (unnormalized).
% L     : truncated order (L <= N).
%
% Output:
% F     : Inverse Fourier-Legendre transform [F1 ... FM] at M angles.
% H     : matrix of the inverse Fourier-Legendre transform. Size(H) = M x N.
%
% See also flt, sht, isht.

% César D. Salvador
% salvador@perception3d.com
% https://cesardsalvador.github.io/
% https://www.perception3d.com/
% February 11, 2024


function [F, varargout] = iflt(Fn, theta, varargin)
[~, Nb] = size(Fn); % Nf : length of frequency domain, Nb : number of coeff.
N = Nb;
F = 0;
H = zeros(length(theta), N);
if (nargin < 3)
    for n = 0:N-1
            hn = pnm(n, 0, sin(theta)).';
            F = F + Fn(:, n+1)*hn;
            H(:, n+1) = hn;
    end
elseif (nargin == 3) && ( strcmp(varargin{1}, 'sch') || ...
                          strcmp(varargin{1}, 'norm') || ...       
                          strcmp(varargin{1}, 'none') )
    for n = 0:N-1
            hn = pnm(n, 0, sin(theta), varargin{1}).';
            F = F + Fn(:, n+1)*hn;
            H(:, n+1) = hn;
    end
elseif (nargin == 4) && ( strcmp(varargin{1}, 'sch') || ...
                          strcmp(varargin{1}, 'norm') || ...       
                          strcmp(varargin{1}, 'none') ) ...
                     && (varargin{2} > 0)
    N = min(Nb, varargin{2});
    for n = 0:N-1
            hn = pnm(n, 0, sin(theta), varargin{1}).';
            F = F + Fn(:, n+1)*hn;
            H(:, n+1) = hn;
    end
end
varargout(1) = {H};

