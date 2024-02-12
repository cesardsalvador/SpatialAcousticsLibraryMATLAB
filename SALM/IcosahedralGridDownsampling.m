%% Indices for downsampling of icosahedral grids based on dyadic edge subdivision

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

close all
clear all
clc

%% load highest resolution icosahedral grid (xh)
ico_edge_div_high = 32;                                 % subdivision factor of icosahedron's edges
number_of_points_high = 10*ico_edge_div_high^2+2;       % number of points
file_name = ['icolr', num2str(ico_edge_div_high), '.mat'];
load(file_name);
xh = x;

%% load low resolution icosahedral grid (xl)
ico_edge_div_low = 1;                                  % subdivision factor of icosahedron's edges
number_of_points_low = 10*ico_edge_div_low^2+2;         % number of points
file_name = ['icolr', num2str(ico_edge_div_low), '.mat'];
load(file_name);
xl = x;

%% Indices for downsampling from xh to xl
[~, indico] = max(xh*xl');
error = max(abs(xh(indico, :) - xl))
% index_file_name = ['indicolr', num2str(ico_edge_div_high), 'to', num2str(ico_edge_div_low), '.mat'];
% save(index_file_name, 'indico');
