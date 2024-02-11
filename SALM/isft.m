% Inverse spherical Fourier Transform.
%    F = isft(Fnm, x)
%    F = isft(Fnm, x, value)
%    F = isft(Fnm, x, value, norm)
%    F = isft(Fnm, x, value, norm, L)
%    [F, Y] = isft(Fnm, x, ...)
%
% Input:
%  Fnm   : Spherical spectrum [F00, ..., FNN] of order N and length (N+1)^2.
%  x     : Cartesian positions x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)]
%  value : 'real' (default real-valued output) or 'complex' (complex-valued output)
%  norm  : 'norm' (default full norm), 'sch' (Schmidt seminorm), or 'unnorm' (unnormalized).
%  L     : Truncated order (L <= N).
%
% Output:
%  F     : Inverse spherical Fourier transform [F1, ..., FP] at P positions.
%  Y     : Matrix of the inverse spherical Fourier transform. Size(H) = P x (N+1)^2.
%
% See also sft, ynm, pnm, bmf, dvf.

% Cesar D. Salvador
% October 14, 2020
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/

% Reference
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.

function [F, varargout] = isft(Fnm, x, varargin)
Q = size(Fnm, 2);       % Q: Number of SFT coefficients
N = ceil(sqrt(Q)-1);    % N: Order of the SFT
F = 0;
Y = zeros(size(x, 1), (N+1)^2);
if nargin < 3
    for n = 0:N
        for m = -n:n
            y = ynm(n, m, x)';
            F = F + Fnm(:, n^2+n+m+1)*y;
            Y(:, n^2+n+m+1) = y;
        end
    end
elseif (nargin == 3) ...
        && any(strcmp(varargin{1}, {'real', 'complex'}))
    for n = 0:N
        for m = -n:n
            y = ynm(n, m, x, varargin{1}).';
            F = F + Fnm(:, n^2+n+m+1)*y;
            Y(:, n^2+n+m+1) = y;
        end
    end
elseif (nargin == 4) ...
        && any(strcmp(varargin{1}, {'real', 'complex'})) ...
        && any(strcmp(varargin{2}, {'norm', 'sch', 'unnorm'}))
    for n = 0:N
        for m = -n:n
            y = ynm(n, m, x, varargin{1}, varargin{2}).';
            F = F + Fnm(:, n^2+n+m+1)*y;
            Y(:, n^2+n+m+1) = y;
        end
    end
elseif (nargin == 5) ...
        && any(strcmp(varargin{1}, {'real', 'complex'})) ...
        && any(strcmp(varargin{2}, {'norm', 'sch', 'unnorm'})) ...
        && (varargin{3} > 0)
    for n = 0:min(N, varargin{3})
        for m = -n:n
            y = ynm(n, m, x, varargin{1}, varargin{2}).';
            F = F + Fnm(:, n^2+n+m+1)*y;
            Y(:, n^2+n+m+1) = y;
        end
    end
else
    disp('Output zero. Value must be REAL or COMPLEX. Norm must be NORM, SCH or UNNORM.')
end
varargout(1) = {Y};
