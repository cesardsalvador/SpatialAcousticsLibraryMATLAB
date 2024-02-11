% Maximum order of a spherical Fourier transform
%    N = sftmaxorder(a, b, k, error)
%
% Input:
% a        : radius of the measurement surface (m)
% b        : radius of the sound source (m)
% k        : wave number (2*pi*f/c)
% err      : specified order truncation error (e.g. 1e-3)
%
% Output:
% N        : maximum order.
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


function N = sftmaxorder(a, b, k, err)
ka = k*a;
r = b/a;
p = 4;
Ncurve = ka + 0.5 * ( 3*log(1./err) + 0.5*log(ka) ).^(2/3) .* (ka).^(1/3);
Noffset = 1./log(r) .* log( r^(3/2) ./ ( (r-1).^(3/2) * err ) ) + 1;
N = floor( abs( Ncurve.^p + Noffset.^p ).^(1./p) );
