% Rigid sphere boundary filtering for spherical harmonics encoding
%    hn = boundfiltn(n, a, r, rs, f, c, field, scale)
% Input
% n     : order of the decomposition
% a     : radius of the rigid sphere (scalar)
% r     : radius of the measurement sphere r >= a (scalar)
% rs    : radius of the source's position (scalar; used if field='near')
% f     : frequency domain (vector)
% c     : velocity of sound in air
% field : 'far' (farfield source), 'nea' (nearfield source), or 'mod'
% scale : 1 (scaled pressure) or any value (not scaled pressure)
% Output
% hn    : sphere boundary weighting coefficients.

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

% Sphere boundary filter for spherical nearfield acoustic holography
function hn = boundfiltn(n, a, r, rs, f, c, field, scale)
Nf = length(f);
k = 2*pi*f/c;
kind = 1;
jnr = besseljsph(n, k*r);
hnr = besselhsph(n, kind, k*r);
djna = besseljsph(n-1, k*a) - ...
    (n+1)./(k*a).*besseljsph(n, k*a);
dhna = besselhsph(n-1, kind, k.*a) - ...
    (n+1)./(k*a).*besselhsph(n, kind, k*a);
bn = jnr - djna./dhna.*hnr;
if field == 'far'
    bn = bn;
    if scale == 1
        bn = (N+1)^2/(4*pi)*bn;
    end
elseif field == 'nea'
    hnrs = besselhsph(n, kind, k*rs);
%    bn = 4*pi*1j*hnrs.*bn;
    bn = -hnrs.*bn;
    if scale == 1
        bn = (N+1)^2*k*rs./(4*pi*exp(1j*k*rs)).*bn;
    end
else
    disp('field must be far or near')
    return
end
hn = 1./bn;
return

% Spherical Bessel function
function h = besseljsph(m, x)
h = sqrt(pi./(2*x)).*besselj(m+0.5, x);

% Spherical Hankel function
function h = besselhsph(m, kind, x)
h = sqrt(pi./(2*x)).*besselh(m+0.5, kind, x);

