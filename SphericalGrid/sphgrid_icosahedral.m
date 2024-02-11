% Based the parametrised icosahedron showing some cube symmetries
% fon GeoGebra
% https://www.geogebra.org/m/qRm7fGdS

close all
clear all
clc

cmap = flipud(colormap(gray(256)));
font_size = 10;
set(0, 'defaultfigurecolor', 'w', ...
    'defaultfigurecolormap', cmap, ...
    'defaultaxesfontname',  'times', ...
    'defaulttextfontSize', font_size, ...
    'defaultaxesfontsize', font_size, ...
    'defaultaxeslayer', 'top')

r = 85;
divline = 2;
number_of_directions = (divline^2)*10+2;

a = 1 / sqrt(5);
b = (1 - a) / 2;
c = (1 + a) / 2;
d = sqrt(b);
e = sqrt(c);

ico = [0 0 1;
        0 2*a a;
        e b a;
        d -c a;
        -d -c a;
        -e b a;
        d c -a;
        e -b -a;
        0 -2*a -a;
        -e -b -a;
        -d c -a;
        0 0 -1];

%thetax = acosd(sqrt(5)/3) - 10; thetay = 0; thetaz = -18*5;
thetax = 0; thetay = 0; thetaz = -18*5;
Rx = [1 0 0; 0 cosd(thetax) -sind(thetax); 0 sind(thetax) cosd(thetax)];
Ry = [cosd(thetay) 0 sind(thetay); 0 1 0; -sind(thetay) 0 cosd(thetay)];
Rz = [cosd(thetaz) -sind(thetaz) 0; sind(thetaz) cosd(thetaz) 0; 0 0 1];
%ico = (Rz*Rx*ico')';

poly = ico .* r;


neighbor_length = poly(2,:) - poly(1,:);
neighbor_length = sqrt(neighbor_length(1)^2 + neighbor_length(2)^2 + neighbor_length(3)^2);

error = neighbor_length*0.2;


for i=1:length(poly)-1
    for j=i+1:length(poly)
        clear judge
        judge = poly(j,:) - poly(i,:);
        judge = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);
        if abs(judge - neighbor_length) < error
            plot3(poly([i j], 1), poly([i j], 2), poly([i j], 3));
            hold on
            neighbor(i, j) = 1;
        end
    end
end

%view([-1 0 0]);
%axis([-100 100 -100 100 -100 100]);
zoom(0.7), grid on, box on, axis equal, axis vis3d, xlabel('x'), ylabel('y'), zlabel('z')
title(['20 polyhedron']);
view([45 45])
%axis square;


clear neighbor_length



k=1;


for i=1:length(poly)-1
    
    poly2(k,:) = poly(i,:);
    k = k+1;
    for j=i+1:length(poly)
        
        if neighbor(i,j) == 1
            
            new_pt_len = (poly(j, :) - poly(i, :));
            
            new_pt_len_abs = sqrt(new_pt_len(1)^2 + new_pt_len(2)^2 + new_pt_len(3)^2);

            a = new_pt_len_abs;
            
            A = acos((r^2+r^2-a^2)/(2*r*r));
            BC = (pi-A)/2;
            
            for l = 1:divline-1
                a_len = r*sin(A*l/divline)/sin(pi-BC-A*l/divline);
                new_pt = poly(i, :) + new_pt_len*(a_len/new_pt_len_abs);
            
                poly2(k,:) = new_pt;
                k = k+1;
            end
        end
    end
end

poly2(k,:) = poly(length(poly),:);


for i=1:length(poly2)
    [THETA PHI R] = cart2sph(poly2(i,1),poly2(i,2),poly2(i,3));
    poly2_sph(i,:) = [THETA poly2(i,3)];
end


[Y I] = sortrows(poly2_sph,[-2, 1]);

for i=1:length(poly2_sph)
    test(i,:) = poly2(I(i),:);
end
value = test;


clear poly poly2


k=1;
elev=1;
count=2;
median_count=0;
for i=1:length(value)
    poly(k,:) = value(i,:);
    k=k+1;
    if i~=length(value)
        if elev < divline
            judge = value(i+1,:) - value(i,:);
            judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);
            if abs(value(i,3) - value(i+1,3)) < 0.1
                a = judge2;
                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;

                for l = 1:elev-1
                    a_len = r*sin(A*l/elev)/sin(pi-BC-A*l/elev);
                    new_pt = value(i, :) + judge*(a_len/judge2);

                    poly(k,:) = new_pt;
                    k = k+1;
                end
                
            elseif abs(value(i,3) - value(i+1,3)) > 1 && i ~=1
                judge = value(count,:) - value(i,:);
                judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);
                a = judge2;
                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;
                for l = 1:elev-1
                    a_len = r*sin(A*l/elev)/sin(pi-BC-A*l/elev);
                    new_pt = value(i, :) + judge*(a_len/judge2);

                    poly(k,:) = new_pt;
                    k = k+1;
                end
                count=i+1;
                elev = elev + 1;
            end

        elseif elev == divline
            if abs(value(i,3) - value(i+1,3)) > 1 && i ~=1
                count=i+1;
                elev = elev + 1;
            end

        elseif elev < 2*divline
            median_count = median_count+1;
            judge = value(i+1,:) - value(i,:);
            judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);

            if abs(value(i,3) - value(i+1,3)) < 0.1
                a = judge2;
                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;
                if rem(floor(2*(elev-divline)/(divline+1))+median_count, 2)
                    for l = 1:elev-divline-1
                        a_len = r*sin(A*l/(elev-divline))/sin(pi-BC-A*l/(elev-divline));
                        new_pt = value(i, :) + judge*(a_len/judge2);

                        poly(k,:) = new_pt;
                        k = k+1;
                    end
                else                    
                    for l = 1:2*divline-elev-1
                        a_len = r*sin(A*l/(2*divline-elev))/sin(pi-BC-A*l/(2*divline-elev));
                        new_pt = value(i, :) + judge*(a_len/judge2);

                        poly(k,:) = new_pt;
                        k = k+1;
                    end
                end

            elseif abs(value(i,3) - value(i+1,3)) > 1 && i ~=1
                judge = value(count,:) - value(i,:);
                judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);
                a = judge2;
                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;
                if rem(floor(2*(elev-divline)/(divline+1))+median_count, 2)
                    for l = 1:elev-divline-1
                        a_len = r*sin(A*l/(elev-divline))/sin(pi-BC-A*l/(elev-divline));
                        new_pt = value(i, :) + judge*(a_len/judge2);

                        poly(k,:) = new_pt;
                        k = k+1;
                    end
                else                    
                    for l = 1:2*divline-elev-1
                        a_len = r*sin(A*l/(2*divline-elev))/sin(pi-BC-A*l/(2*divline-elev));
                        new_pt = value(i, :) + judge*(a_len/judge2);

                        poly(k,:) = new_pt;
                        k = k+1;
                    end
                end
                count=i+1;
                elev = elev + 1;
                median_count=0;
            end
            
        elseif elev == 2*divline
            if abs(value(i,3) - value(i+1,3)) > 1 && i ~=1
                count=i+1;
                elev = elev + 1;
            end

        elseif elev < 3*divline

            judge = value(i+1,:) - value(i,:);
            judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);

            if abs(value(i,3) - value(i+1,3)) < 0.1

                a = judge2;

                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;

                for l = 1:(3*divline - elev) -1

                    a_len = r*sin(A*l/(3*divline - elev))/sin(pi-BC-A*l/(3*divline - elev));
                    new_pt = value(i, :) + judge*(a_len/judge2);

                    poly(k,:) = new_pt;
                    k = k+1;
                end

            elseif abs(value(i,3) - value(i+1,3)) > 1 && i ~=1
                judge = value(count,:) - value(i,:);
                judge2 = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);
                a = judge2;
                A = acos((r^2+r^2-a^2)/(2*r*r));
                BC = (pi-A)/2;
                for l = 1:(3*divline - elev)-1
                    a_len = r*sin(A*l/(3*divline - elev))/sin(pi-BC-A*l/(3*divline - elev));
                    new_pt = value(i, :) + judge*(a_len/judge2);

                    poly(k,:) = new_pt;
                    k = k+1;
                end
                count=i+1;
                elev = elev + 1;
            end

            
        end
    end
end


k=k-1;



for i=1:k
    [THETA PHI R] = cart2sph(poly(i,1),poly(i,2),poly(i,3));

    [X Y Z] = sph2cart(THETA, PHI, r);

    poly2(i,:) = [X Y Z];
end


neighbor_length = poly2(2,:) - poly2(1,:);
neighbor_length = sqrt(neighbor_length(1)^2 + neighbor_length(2)^2 + neighbor_length(3)^2);

error = neighbor_length*0.4;
side_count = 1;

for i=1:length(poly2)-1
    for j=i+1:length(poly2)
        clear judge

        judge = poly2(j,:) - poly2(i,:);
        judge = sqrt(judge(1)^2 + judge(2)^2 + judge(3)^2);

        if abs(judge - neighbor_length) < error
            neighbor2(i, j) = 1;
            side_length(side_count) = judge;
            side_length_error(side_count) = abs(judge - neighbor_length);
            side_count = side_count + 1;
        end
    end
end


clear value

[s  I] = sort(poly2(:,3),'ascend');

for i=1:length(I)
    sortpoly(i,:) = poly2(I(i),:);
end
poly2 = sortpoly;

for i=1:length(sortpoly)
    [THETA PHI R] = cart2sph(sortpoly(i,1), sortpoly(i,2), sortpoly(i,3));
    value(length(sortpoly)+1-i,:) = [THETA PHI];
end

dist = 1;
[x(:, 1), x(:, 2), x(:, 3)] = sph2cart(value(:, 1), value(:, 2), dist);
weight = quadweight(x)'/(4*pi);
msh = convhull(x);

% Spherical harmonic analysis
thres = 1e-4;                   % minimum threshold for plot scale
thres_log = log10(thres);       % threshold in log-scale
ctick_log = thres_log:0;
ctick = 10.^ctick_log;
val = 'complex';                % values of spherical harmonics (real or complex)
nrm = 'norm';                   % normalization type of spherical harmonics
Nx = length(x);
ord = floor(sqrt(Nx) - 1);
if ord > 5, ordinit = 5; else ordinit = 1; end
Y = zeros(Nx, (ord + 1)^2);
for n = 0:ord
    for m=-n:n
        Y(:, n^2+n+m+1) = ynm(n, m, x, val, nrm);
    end
end
D = Y'*diag(weight)*Y;
%D = Y'*Y/Nx;
E = eye(size(D)) - D;

figure
subplot(231)
plot3(x(:, 1), x(:, 2), x(:, 3), 'r.')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis equal, axis tight
view([90 0])

subplot(232)
plot(value(:, 1)*180/pi, value(:, 2)*180/pi, 'r.')
xlabel('Azimuth')
ylabel('Elevation')
set(gca, 'xtick', -180:90:180, 'xticklabel', -180:90:180)
set(gca, 'ytick', -90:45:90, 'yticklabel', -90:45:90)
grid on
axis equal, axis tight

subplot(233)
V = log10(abs(D)); V(V < thres_log) = nan;
imagesc(0:size(D, 1)-1, 0:size(D, 2)-1, V)
axis xy
set(gca, 'xtick', (ordinit:ord).^2+(ordinit:ord), 'xticklabel', ordinit:ord)
set(gca, 'ytick', (ordinit:ord).^2+(ordinit:ord), 'yticklabel', ordinit:ord)
colorbar
set(gca, 'clim', [thres_log 0])
title_str= '$\epsilon^{(Q)}_{nmn^{\prime}m^{\prime}}$';
% title(title_str, 'interpreter', 'latex')
xlabel('$n^{\prime}$', 'interpreter', 'latex')
ylabel('$n$', 'interpreter', 'latex')
hcb = colorbar('eastoutside', ...
    'ytick', ctick_log, 'yticklabel', ctick, 'ylim', [min(ctick_log) max(ctick_log)]);
% set(hcb, 'ytickLabel', cellstr(num2str(reshape(get(hcb, 'ytick'), [], 1),'%0.3f')) )
%set(get(hcb, 'xlabel'), 'String', title_str, 'interpreter', 'latex', 'fontsize', font_size)
%colormap(cmap)
axis equal
axis tight
grid off
axis on
box on
set(gcf, 'paperpositionmode', 'auto')
% fig_name = ['fig_orthonormality_mic_', sampling_mic, '_ord', num2str(ord_mic), fig_nameext];
% print(fig_format, '-tiff', fig_res, [fig_folder, fig_name])


x = (Rz*Rx*x')';
[value(:, 1), value(:, 2), ~] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
weight = quadweight(x)'/(4*pi);
sum(weight)
msh = convhull(x);

% Spherical harmonic analysis
Nx = length(x);
ord = floor(sqrt(Nx) - 1);
if ord > 5, ordinit = 5; else ordinit = 1; end
Y = zeros(Nx, (ord + 1)^2);
for n = 0:ord
    for m=-n:n
        Y(:, n^2+n+m+1) = ynm(n, m, x, val, nrm);
    end
end
D = Y'*diag(weight)*Y;
%D = Y'*Y/Nx;
E = eye(size(D)) - D;

subplot(234)
plot3(x(:, 1), x(:, 2), x(:, 3), 'r.')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis equal, axis tight
view([90 0])

subplot(235)
plot(value(:, 1)*180/pi, value(:, 2)*180/pi, 'r.')
xlabel('Azimuth')
ylabel('Elevation')
set(gca, 'xtick', -180:90:180, 'xticklabel', -180:90:180)
set(gca, 'ytick', -90:45:90, 'yticklabel', -90:45:90)
grid on
axis equal, axis tight

subplot(236)
V = log10(abs(D)); V(V < thres_log) = nan;
imagesc(0:size(D, 1)-1, 0:size(D, 2)-1, V)
axis xy
set(gca, 'xtick', (ordinit:ord).^2+(ordinit:ord), 'xticklabel', ordinit:ord)
set(gca, 'ytick', (ordinit:ord).^2+(ordinit:ord), 'yticklabel', ordinit:ord)
colorbar
set(gca, 'clim', [thres_log 0])
title_str= '$\epsilon^{(Q)}_{nmn^{\prime}m^{\prime}}$';
% title(title_str, 'interpreter', 'latex')
xlabel('$n^{\prime}$', 'interpreter', 'latex')
ylabel('$n$', 'interpreter', 'latex')
hcb = colorbar('eastoutside', ...
    'ytick', ctick_log, 'yticklabel', ctick, 'ylim', [min(ctick_log) max(ctick_log)]);
% set(hcb, 'ytickLabel', cellstr(num2str(reshape(get(hcb, 'ytick'), [], 1),'%0.3f')) )
%set(get(hcb, 'xlabel'), 'String', title_str, 'interpreter', 'latex', 'fontsize', font_size)
%colormap(cmap)
axis equal
axis tight
grid off
axis on
box on
set(gcf, 'paperpositionmode', 'auto')
% fig_name = ['fig_orthonormality_mic_', sampling_mic, '_ord', num2str(ord_mic), fig_nameext];
% print(fig_format, '-tiff', fig_res, [fig_folder, fig_name])


% file_name = ['icolr', num2str(divline), '.mat'];
% save(['D:\prog\matlab\data\ico\hrir\', file_name], 'x', 'msh', 'weight');

disp(['number of points: ' num2str(number_of_directions)]);
a = x * x';
a(abs(a) >= 1) = 0;
disp(['maximum angle between adjacent points: ' num2str(max(min(acosd(a))))]);
