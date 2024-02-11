%% Indices for downsampling of icosahedral grids based on dyadic edge subdivision

% César D. Salvador
% cesardsalvador@gmail.com
% http://www.ais.riec.tohoku.ac.jp/~salvador/

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
