% Spherical harmonic function of order n and degree m
% y = ynm(n, m, x)
% y = ynm(n, m, x, value)
% y = ynm(n, m, x, value, norm)
%
% Input:
%  n     : Order of the spherical harmonic, n = 0, 1, 2, ...
%  m     : Degree of the spherical harmonic, abs(m) <= n.
%  x     : Cartesian positions x = [x(1:P, 1) x(1:P, 2) x(1:P, 3)] on a spherical surface
%  value : 'real' (default real-valued output) or 'complex' (complex-valued output)
%  norm  : 'norm' (default full norm), 'sch' (Schmidt seminorm), or 'unnorm' (unnormalized).
%
% Output:
%  y     : Spherical harmonic function of order n and degree m. Size P x 1.
%
% See also pnm, sft, isft.

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


function y = ynm(n, m, x, varargin)
[azim, elev, ~] = cart2sph(x(:,1), x(:,2), x(:,3));
absm = abs(m);
if nargin < 4
    if m == 0
        fAzim = 1;
    elseif m > 0
        fAzim = sqrt(2) * cos(m*azim);
    else
        fAzim = sqrt(2) * sin(absm*azim);
    end
    Nnm = (-1)^absm * sqrt((2*n+1)/(4*pi)*factorial(n-absm)/factorial(n+absm));
elseif (nargin == 4) ...
        && any(strcmp(varargin{1}, {'real', 'complex'}))
    if strcmp(varargin{1}, 'real')
        if m == 0
            fAzim = 1;
        elseif m > 0
            fAzim = sqrt(2) * cos(m*azim);
        else
            fAzim = sqrt(2) * sin(absm*azim);
        end
    elseif strcmp(varargin{1}, 'complex')
        fAzim = exp(1j*m*azim);
    else
        fAzim = 0;
        disp('Output zero. Value must be REAL or COMPLEX.')
    end
    Nnm = (-1)^absm * sqrt((2*n+1)/(4*pi)*factorial(n-absm)/factorial(n+absm));
elseif (nargin == 5) ...
        && any(strcmp(varargin{1}, {'real', 'complex'})) ...
        && any(strcmp(varargin{2}, {'norm', 'sch', 'unnorm'}))
    switch varargin{1}
        case 'real'
            if m == 0
                fAzim = 1;
            elseif m > 0
                fAzim = sqrt(2) * cos(m*azim);
            else
                fAzim = sqrt(2) * sin(absm*azim);
            end
        case 'complex'
            fAzim = exp(1j*m*azim);
        otherwise
            fAzim = 0;
            disp('Output zero. Value must be REAL or COMPLEX.')
    end
    switch varargin{2}
        case 'norm'
            Nnm = (-1)^absm * sqrt((2*n+1)/(4*pi)*factorial(n-absm)/factorial(n+absm));
        case 'sch'
            if m == 0
                Nnm = 1;
            else
                Nnm = (-1)^absm * sqrt(2*factorial(n-absm)/factorial(n+absm));
            end
        case 'unnorm'
            Nnm = 1;
        otherwise
            Nnm = 0;
            disp('Output zero. Norm must be NORM, SCH or UNNORM.')
    end
else
    fAzim = 0;
    Nnm = 0;
    disp('Output zero. Value must be REAL or COMPLEX. Norm must be NORM, SCH or UNNORM.')
end
fElev = pnm(n, absm, sin(elev));
y = Nnm * fAzim .* fElev;
