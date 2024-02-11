% Transform-domain boundary matching filters (BMFs) for spherical arrays
%        Fn = bmf(a, b, k, n, kind, type)
%        Fn = bmf(a, b, k, n, kind, type, alpha)
%        [Fn, rhon] = bmf(a, b, k, n, kind, 'ROlim', alpha)
%
% INPUT
%	a:     radius of recording boundary
%	b:     radius of reproduction boundary (input any number if type is 'ri')
%	k:     vector of wave numbers
%	n:     order of spherical Hankel functions
%	kind:  kind of spherical Hankel functions (1 or 2)
%   type:  type pf BMFs
%          'RI' (Rigid recording boundary to Infinite reproduction boundary)
%          'RF' (Rigid recording boundary to Finite reproduction boundary)
%          'RO' (Rigid recording boundary to Open reproduction boundary)
%          'RItrunc1' (truncated band-limited BMFs of type 'RI', constant reproduction radius)
%          'RFtrunc1' (truncated band-limited BMFs of type 'RF', constant reproduction radius)
%          'RItrunc2' (truncated band-limited BMFs of type 'RI', constant order truncation error)
%          'RFtrunc2' (truncated band-limited BMFs of type 'RF', constant order truncation error)
%          'RIreg' (regularization-based band-limited BMFs of type 'RI')
%          'RFreg' (regularization-based band-limited BMFs of type 'RF')
%          'ROlim' (limited magnitude band-limited BMFs of type 'RO')
%   alpha: scalar parameter for band-limitation
%          If 'RItrunc1' or 'RFtrunc1', alpha is the reproduction radius (e.g. alpha = b)
%          If 'RItrunc2' or 'RFtrunc2', alpha is the truncation error (e.g. alpha = 5e-2)
%          If 'RIreg' or 'RFreg', alpha is the regularization parameter (e.g. alpha = 5e-3)
%          If 'ROlim', alpha is the maximum magnitude squared in decibell (e.g. alpha = 40)
%
% OUTPUT
%	Fn:    boundary matching filter of order n
%	rhon:  modified reproduction radius for BMFs of type 'RO' and order n.
%
% See also obf, besselhsph, dbesselhsph, sft, isft.

% Cesar D. Salvador
% https://cesardsalvador.github.io
% February 18, 2021

% Reference and citation
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.
%     DOI: 10.1109/TASLP.2017.2778562
% [2] C. D. Salvador et al., “Design theory for binaural synthesis:
%     Combining microphone array recordings and head-related transfer function datasets,”
%     Acoust. Sci. Technol., vol. 38, no. 2, pp. 51–62, Mar. 2017.
%     DOI: 10.1250/ast.38.51

function [Fn, varargout] = bmf(a, b, k, n, kind, type, varargin)
dhna = dbesselhsph(n, kind, k*a);
hnb = besselhsph(n, kind, k*b);
switch type
    case 'RI'
        Fn = (k * a).^2 .* dhna ./ ( 4*pi * (1j)^(n+1) );
    case 'RF'
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb );
    case 'RO'
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb );
        Fn = obf(b, k, n, kind) .* Fn;
    case 'RItrunc1'
        mask = n <= floor( exp(1)/2 * k .* varargin{1} );
        Fn = (k * a).^2 .* dhna ./ ( 4*pi * (1j)^(n+1) ) .* mask;
    case 'RFtrunc1'
        mask = n <= floor( k .* varargin{1} );
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb ) .* mask;
    case 'RItrunc2'
        mask = n <= sftmaxorder(a, b, k, varargin{1});
        Fn = (k * a).^2 .* dhna ./ ( 4*pi * (1j)^(n+1) ) .* mask;
    case 'RFtrunc2'
        mask = n <= sftmaxorder(a, b, k, varargin{1});
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb ) .* mask;
    case 'RIreg'
        Fn = -(k * a).^2 .* dhna ./ ( 4*pi * (1j)^(n+1) );
        Fn = Fn ./ ( 1 + varargin{1}.^2 .* abs(Fn).^2 );
    case 'RFreg'
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb );
        Fn = Fn ./ ( 1 + varargin{1}.^2 .* abs(Fn).^2 );
    case 'ROlim'
        Fn = -(k * a^2) .* dhna ./ ( 4*pi * hnb );
        ind_k_gamma = find( 20*log10(abs(Fn)) <= varargin{1}, 1, 'first' );
        if ind_k_gamma == 1
            rhon = b;
        else
            rhon = 0.98 * (n + 0.5) / k(ind_k_gamma);
            Fn = obf(rhon, k, n, kind) .* Fn;
        end
        varargout{1} = min(rhon, b);
    otherwise
        Fn = 0;
        disp('Fn = 0. Invalid type of BMFs.')
end
