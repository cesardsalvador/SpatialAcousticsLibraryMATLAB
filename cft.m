% Circular Fourier transform (CFT) along a dimension
%    Fm = cft(F, dim)
%
% INPUT
%   F:      Function in circular domain
%   dim:    dimension along which CFT is applied
%
% OUTPUT
%	Fm:     CFT of F with symmetric ordering m = [-M/2:M/2-1]
%
% See also icft, flt, iflt, dvf.

% Cesar D. Salvador
% cesardsalvador@gmail.com

% Reference and citation
% [1] C. D. Salvador et al., “Distance-varying filters to synthesize
%     head-related transfer functions in the horizontal plane
%     from circular boundary values,” Acoust. Sci. Technol.,
%     vol. 38, no. 1, pp. 1–13, Jan. 2017.
%     DOI: 10.1250/ast.38.1

function Fm = cft(F, dim)
Fm = fftshift(fft(fftshift(F, dim), [], dim), dim);
