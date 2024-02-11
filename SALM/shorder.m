% Maximum order of a spherical harmonic decomposition.
%    N = shorder(rm, rs, k, error)
%
% Input:
% rm        : radius of the measurement surface [m].
% rs        : radius of the sound source [m].
% k         : wave number (2*pi*f/c).
% err       : desired approximation error.
%
% Output:
% N         : maximum order.
%
% See also sst, isst, isht, sht.

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


function N = shorder(rm, rs, k, err)
%a = rm/rs;
%Noffset = log(3/2*(1-a)*err) / log(a) + 1;
a = rs/rm;
Noffset = 1/log(a) * log( a^(3/2) / ( (a-1)^(3/2) * err ) ) + 1;
Ncurve = k*rm + 0.5 * ( 3*log(1./err) + 0.5*log(k*rm) ).^(2/3) .* ...
                (k*rm).^(1/3);
p = 4;
N = ceil( abs( Noffset.^p + Ncurve.^p ).^(1./p) );
