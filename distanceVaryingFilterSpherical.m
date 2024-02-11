% Distance-varying filter (DVF) of order n in the spherical Fourier transform (SFT) domain.
% DVF extrapolates from distance a to distance b in the frequency domain.
% DVF requires a spherical set of acoustic data.
%
%  Dn = distanceVaryingFilterSpherical(a, b, k, regParam)
%
% Input:
%  a        : Initial distance in meters.
%  b        : Final distance in meters.
%  k        : Wave number [k(1); k(2); ... ; k(Ns/2+1)], Ns is number of samples in time.
%  regParam : Regularization parameter to avoid excessive gains (e.g., 1e-2).
%
% Output:
%  Dn       : Distance varying filter of order n. Column vector of lenght Ns/2+1.
%
% See also besselhsph, sft, isft

% November 22, 2020
% Cesar D. Salvador
% cesardsalvador@gmail.com
% https://cesardsalvador.github.io/

% Reference
% [1] R. Duraiswami et al., “Interpolation and range extrapolation of HRTFs,”
% in Proc. IEEE ICASSP, May 2004, vol. 4, pp. 45–48, doi: 10.1109/ICASSP.2004.1326759.
% [2] M. Pollow et al., “Calculation of head-related transfer functions for arbitrary field
% points using spherical harmonics,” Acta Acust. United Ac., vol. 98, no. 1, pp. 72–82, Jan. 2012.
% [3] C. D. Salvador et al., “Distance-varying filters to synthesize head-related
% transfer functions in the horizontal plane from circular boundary values,”
% Acoust. Sci. Technol., vol. 38, no. 1, pp. 1–13, Jan. 2017, doi: 10.1250/ast.38.1.

function Dn = distanceVaryingFilterSpherical(a, b, n, k, regParam)
Ms = length(k);                                 % Ms = Ns/2+1, Ns = Number of samples in time.
Dn = zeros(Ms, 1);                              % Initialization of DVF
if b < a
    kind = 1;                                   % First kind of spherical Hankel functions
else
    kind = 2;                                   % Second kind of spherical Hankel functions
end
hna = besselhsph(n, kind, k(2:end)*a);          % Spherical Hankel function at a
hnb = besselhsph(n, kind, k(2:end)*b);          % Spherical Hankel function at b
Dn(2:end) = hnb ./ hna;                         % DVF from a to b
WnReg = 1 ./ (1 + regParam^2*abs(Dn).^2);       % Window for regularization
Dn = WnReg .* Dn;                               % Regularized DVF
