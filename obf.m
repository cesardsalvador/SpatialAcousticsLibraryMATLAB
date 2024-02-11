% Spherical open boundary filters for angular band-limiting
% array
%        Fn = obf(b, k, n, kind)
% Input:
%	b:     radius of open sphere
%	k:     wave number
%	n:     order of spherical Hankel functions
%	kind:  kind of spherical Hankel functions
% Output
%	Fn:    open boundary filter
%
% See also bmf, besselhsph, dbesselhsph, sft, isft.

% Cesar D. Salvador
% cesardsalvador@gmail.com

% Reference and citation
% [1] C. D. Salvador et al., “Boundary matching filters for spherical
%     microphone and loudspeaker arrays,” IEEE/ACM Trans. Audio, Speech, Language Process.,
%     vol. 26, no. 3, pp. 461–474, Mar. 2018.
%     DOI: 10.1109/TASLP.2017.2778562

function Fn = obf(b, k, n, kind)
hnb = besselhsph(n, mod(kind,2)+1, k*b);
dhnb = dbesselhsph(n, kind, k*b);
Fn = abs( 1j ./ ( (k*b).^2 .* hnb .* dhnb ) );
