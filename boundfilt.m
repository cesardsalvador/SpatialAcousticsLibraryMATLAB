% Rigid sphere boundary filtering for spherical harmonics encoding
%    hn = boundfilt(N, a, r, rs, f, c, field, scale)
% Input
% N     : order of the decomposition
% a     : radius of the rigid sphere (scalar)
% r     : radius of the measurement sphere r >= a (scalar)
% rs    : radius of the source's position (scalar; used if field='near')
% f     : frequency domain (vector)
% c     : velocity of sound in air
% field : 'far' (farfield source), 'nea' (nearfield source), or 'mod'
% scale : 1 (scaled pressure) or any value (not scaled pressure)
% Output
% hn    : sphere boundary weighting coefficients [h00 ... hNN]

% References
% [1] E. G. Williams, Fourier Acoustics: Sound Radiation and
% Nearfield Acoustical Holography. London, UK: Academic Press, 1999.
% [2] J. J. Bowman, T. B. A. Senior, and P. L. E. Uslenghi,
% Electromagnetic and acoustic scattering by simple shapes.
% New York, NY, USA: Hemisphere, 1987.

% Cesar D. Salvador
% cdsalv@gmail.com

% Sphere boundary filter for spherical nearfield acoustic holography
function hn = boundfilt(N, a, r, rs, f, c, field, scale)
Nf = length(f);
k = 2*pi*f/c;
kind = 1;
hn = zeros(Nf, N+1);
for n = 0:N
    jnr = besseljsph(n, k*r);
    hnr = besselhsph(n, kind, k*r);
    djna = besseljsph(n-1, k*a) - ...
        (n+1)./(k*a).*besseljsph(n, k*a);
    dhna = besselhsph(n-1, kind, k.*a) - ...
        (n+1)./(k*a).*besselhsph(n, kind, k*a);
    bn = jnr - djna./dhna.*hnr;
    if field == 'far'
        bn = 4*pi*(1j)^n*bn;
        if scale == 1
            bn = (N+1)^2/(4*pi)*bn;
        end
    elseif field == 'nea'
        hnrs = besselhsph(n, kind, k*rs);
        bn = 4*pi*1j*hnrs.*bn;
        if scale == 1
            bn = (N+1)^2*k*rs./(4*pi*exp(1j*k*rs)).*bn;
        end
    else
        disp('field must be far or near')
        return
    end
    for m = -n:n
        hn(:, n^2+n+m+1) = 1./bn;
    end
end
return

% Spherical Bessel function
function h = besseljsph(m, x)
h = sqrt(pi./(2*x)).*besselj(m+0.5, x);

% Spherical Hankel function
function h = besselhsph(m, kind, x)
h = sqrt(pi./(2*x)).*besselh(m+0.5, kind, x);

