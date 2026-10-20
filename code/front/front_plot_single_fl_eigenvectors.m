function front_plot_single_fl_eigenvectors(data, front, outputFile, ...
    globalMaxValue, maxHalfLengthKm, rootDegree, rankValue)
%FRONT_PLOT_SINGLE_FL_EIGENVECTORS Plot one FL and its eigendirections.

fig = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[1 1 19 18]);
cleanupObj = onCleanup(@() close(fig));
ax = axes(fig);
imagesc(ax,data.lonir(1,:),data.latir(:,1),data.M2);
set(ax,'YDir','normal','SortMethod','childorder'); hold(ax,'on');
colormap(ax,summer(256)); cb = colorbar(ax); cb.Label.String = 'M2';

ridge = front.ridge;
hFL = plot(ax,ridge.x,ridge.y,'k-','LineWidth',2.5,'DisplayName','FL');
colors = {[0.68 0.68 0.68],[0.42 0.42 0.42], ...
          [0 0.25 0.95],[0.95 0.02 0.02]};
widths = [1.8 1.8 0.9 1.2];
samples = front.eigen_samples;
for j = 1:numel(samples)
    s = samples(j);
    if any(~isfinite(s.AT_values)), continue; end
    plot(ax,s.x,s.y,'ko','MarkerFaceColor','w','MarkerSize',3.0, ...
        'HandleVisibility','off');
    vectors = {s.S_vectors(:,1),s.S_vectors(:,2), ...
        s.AT_vectors(:,1),s.AT_vectors(:,2)};
    values = [s.S_values(1),s.S_values(2),s.AT_values(1),s.AT_values(2)];
    for q = 1:4
        if ~isfinite(values(q)) || abs(values(q)) <= eps, continue; end
        halfLength = maxHalfLengthKm * ...
            (abs(values(q))/globalMaxValue)^(1/rootDegree);
        vector = halfLength .* vectors{q};
        draw_double_arrow(ax,s.x,s.y,vector(1),vector(2), ...
            colors{q},widths(q),values(q)<0);
    end
end

marginKm = maxHalfLengthKm + 3;
axis(ax,'equal');
xlim(ax,[max(min(data.lonir(:)),min(ridge.x)-marginKm), ...
         min(max(data.lonir(:)),max(ridge.x)+marginKm)]);
ylim(ax,[max(min(data.latir(:)),min(ridge.y)-marginKm), ...
         min(max(data.latir(:)),max(ridge.y)+marginKm)]);
grid(ax,'on');
set(ax,'FontSize',11,'Box','on','Layer','top');
xlabel(ax,'x (km)'); ylabel(ax,'y (km)');
title(ax,sprintf('Rank %d: FL%d flow field, length %.2f km', ...
    rankValue,front.fl_id,ridge.length_km));
h1 = plot(ax,NaN,NaN,'-','Color',colors{1},'LineWidth',widths(1));
h2 = plot(ax,NaN,NaN,'-','Color',colors{2},'LineWidth',widths(2));
h3 = plot(ax,NaN,NaN,'-','Color',colors{3},'LineWidth',widths(3));
h4 = plot(ax,NaN,NaN,'-','Color',colors{4},'LineWidth',widths(4));
legend(ax,[hFL h1 h2 h3 h4],{'FL','S min','S max','A^T negative','A^T positive'}, ...
    'Location','southoutside','NumColumns',2);
exportgraphics(fig,outputFile,'Resolution',600,'BackgroundColor','white');
end

function draw_double_arrow(ax,x0,y0,ux,uy,colorValue,lineWidth,inward)
if hypot(ux,uy) <= eps, return; end
if inward
    quiver(ax,x0+ux,y0+uy,-ux,-uy,0,'Color',colorValue, ...
        'LineWidth',lineWidth,'MaxHeadSize',0.5,'AutoScale','off', ...
        'HandleVisibility','off');
    quiver(ax,x0-ux,y0-uy,ux,uy,0,'Color',colorValue, ...
        'LineWidth',lineWidth,'MaxHeadSize',0.5,'AutoScale','off', ...
        'HandleVisibility','off');
else
    quiver(ax,x0,y0,ux,uy,0,'Color',colorValue,'LineWidth',lineWidth, ...
        'MaxHeadSize',0.5,'AutoScale','off','HandleVisibility','off');
    quiver(ax,x0,y0,-ux,-uy,0,'Color',colorValue,'LineWidth',lineWidth, ...
        'MaxHeadSize',0.5,'AutoScale','off','HandleVisibility','off');
end
end

