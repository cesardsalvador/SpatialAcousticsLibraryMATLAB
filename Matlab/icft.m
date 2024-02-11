% Inverse circular Fourier transform (ICFT) along a dimension
%    F = icft(Fm, dim)
%    F = icft(Fm, dim, n)
%
% INPUT
%   Fm:     Function in circular Fourier domain with symmetric ordering m = [-M/2:M/2-1]
%   dim:    dimension along which ICFT is applied
%
% OUTPUT
%	F:      ICFT of Fm
%
% See also cft, flt, iflt, dvf.

% Cesar D. Salvador
% cesardsalvador@gmail.com

% Reference and citation
% [1] C. D. Salvador et al., “Distance-varying filters to synthesize
%     head-related transfer functions in the horizontal plane
%     from circular boundary values,” Acoust. Sci. Technol.,
%     vol. 38, no. 1, pp. 1–13, Jan. 2017.
%     DOI: 10.1250/ast.38.1

function F = icft(Fm, dim)
F = fftshift(ifft(fftshift(Fm, dim), [], dim, 'nonsymmetric'), dim);
