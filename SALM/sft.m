% Spherical Fourier transform (SFT).
%    Fnm = sft(F, N, x)
%    Fnm = sft(F, N, x, value)
%    Fnm = sft(F, N, x, value, norm)
%    Fnm = sft(F, N, x, value, norm, pinv, tol)
%    Fnm = sft(F, N, x, value, norm, pinv, tol, quad)
%    [Fnm, G] = sft(F, N, x, ...)
%
% Input:
%  F     : Function F = [F1, ..., FP] evaluated at P positions.
%  N     : Order of the SFT (natural number).
%  x     : Cartesian positions x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)]
%  value : 'real' (default real-valued output) or 'complex' (complex-valued output)
%  norm  : 'norm' (default full norm), 'sch' (Schmidt seminorm), or 'unnorm' (unnormalized).
%  pinv  : 'moore' for Moore-Penrose pseudo-inverse (default), or 'reg' for regularization.
%  tol   : Tolerance or regularization parameter to compute the pseudo-inverse. Positive scalar.
%  quad  : Quadrature weights for integration on the sphere. Vector of size Px1. Sum(quad)=1.
%
% Output:
%  Fnm   : Spherical Fourier transform [F00, ..., FNN] of order N and length (N+1)^2.
%  G     : Matrix of the SFT. Size(G) = (N+1)^2 x P.
%
% See also isft, ynm, pnm, bmf, dvf.

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


function [Fnm, varargout] = sft(F, N, x, varargin)

% Matrix of spherical harmonics (Y)
Y = zeros(size(x, 1), (N+1)^2);
if nargin < 4
    for n = 0:N
        for m = -n:n
            Y(:, n^2+n+m+1) = ynm(n, m, x);
        end
    end
elseif (nargin == 4) ...
        && any(strcmp(varargin{1}, {'real', 'complex'}))
    for n = 0:N
        for m = -n:n
            Y(:, n^2+n+m+1) = ynm(n, m, x, varargin{1});
        end
    end
elseif (nargin > 4) ...
        && any(strcmp(varargin{1}, {'real', 'complex'})) ...
        && any(strcmp(varargin{2}, {'norm', 'sch', 'unnorm'}))
    for n = 0:N
        for m = -n:n
            Y(:, n^2+n+m+1) = ynm(n, m, x, varargin{1}, varargin{2});
        end
    end
else
    disp('Output zero. Value must be REAL or COMPLEX. Norm must be NORM, SCH or UNNORM.')
    Fnm = 0;
    varargout(1) = {[]};
    return
end

% Pseudo-inverse of matrix of spherical harmonics, G = pinv(Y)
if nargin < 6
    G = pinv(Y);
elseif (nargin > 5) && (nargin < 8) && (varargin{4} >= 0) && strcmp(varargin{3}, 'moore')
    G = pinv(Y, varargin{4});
elseif (nargin > 5) && (nargin < 8) && (varargin{4} >= 0) && strcmp(varargin{3}, 'reg')
    [U, S, V] = svd(Y, 'econ');
    s = diag(S);
    Sreginv = diag( abs(s).^2 ./ (abs(s).^2+varargin{4}^2) ./ s );
    G = V * Sreginv * U';
elseif (nargin == 8) && ((sum(varargin{5})-1) < 1e-6) && (varargin{4} >= 0) && strcmp(varargin{3}, 'reg')    
    G = (Y'*diag(varargin{5})*Y + varargin{4}^2*eye((N+1)^2)) \ (Y'*diag(varargin{5}));
else
    disp('Output zero. Tol must be positive, Pinv must be MOORE or REG. Veryfy that sum(quad)==1.')
    G = zeros((N+1)^2, size(x, 1));
end

% Spherical Fourier transform
Fnm = (G*F.').';
varargout(1) = {G};
