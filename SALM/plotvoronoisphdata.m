function plotvoronoisphdata(x, v)
radius = max(abs(x(:)));
[~,~,vCells]=voronoisphere(x');
axes('NextPlot','add');
for j=1:length(x)
    fill3(radius*vCells{j}(1,:), ...
        radius*vCells{j}(2,:), ...
        radius*vCells{j}(3,:), ...
        v(j), ...
        'edgecolor', 'none');
end
axis equal
axis vis3d
grid on
colorbar
zoom(0.5)
