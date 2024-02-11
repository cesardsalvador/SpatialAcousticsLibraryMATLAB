% Associated Legendre polynomial of degree n and order m
% p = pnm(n, m, x)
% p = pnm(n, m, x, nrm)
%
% Inputs:
%  n     : Degree of the associated Legendre polynomial, n = 0, 1, 2, ...
%  m     : Order of the associated Legendre polynomial, abs(m) <= n.
%  x     : Argument (scalar, vector or array).
%  norm  : 'unnorm' (default non-normalized), 'norm' (full norm), or 'sch' (Schmidt seminorm).
%
% Outputs:
%  p     : Associated Legendre polynomial of degree n and order m
%
% See also ynm, sft, isft.

% Cesar D. Salvador
% October 14, 2020
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/

% Reference
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.

function p = pnm(n, m, x, varargin)
absn = abs(n);
absm = abs(m);
if nargin < 4
    if n == 0
        p = legendre(n, x(:));
    elseif n > 0
        pos = legendre(n, x(:));
        if m >= 0
            p = pos(m+1,:);
        else
            p = (-1)^absm*factorial(n-absm)/factorial(n+absm) * pos(absm+1,:);
        end
    else
        neg = legendre(absn-1, x(:));
        if m >= 0
            p = neg(m+1,:);
        else
            p = (-1)^absm*factorial(absn-absm)/factorial(absn+absm) * neg(absm+1,:);
        end
    end
elseif (nargin == 4) && any(strcmp(varargin{1}, {'norm', 'sch', 'unnorm'}))
    if n == 0
        p = legendre(n, x(:), varargin{1});
    elseif n > 0
        pos = legendre(n, x(:), varargin{1});
        if m >= 0
            p = pos(m+1,:);
        else
            p = (-1)^absm*factorial(n-absm)/factorial(n+absm) * pos(absm+1,:);
        end
    else
        neg = legendre(absn-1, x(:), varargin{1});
        if m >= 0
            p = neg(m+1,:);
        else
            p = (-1)^absm*factorial(absn-absm)/factorial(absn+absm) * neg(absm+1,:);
        end
    end
else
    disp('Output zero. Norm must be UNNORM, NORM, or SCH')
    p = zeros(size(x));
end
p = reshape(p, size(x));
