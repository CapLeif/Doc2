function front_plot_curated_single_fl(data,front,outputFile,globalMax,maxHalfLengthKm,rootDegree)
%FRONT_PLOT_CURATED_SINGLE_FL Plot one filtered FL with matching arrows.

fig=figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 20 18]);
cleanupObj=onCleanup(@() close(fig)); ax=axes(fig);
imagesc(ax,data.lonir(1,:),data.latir(:,1),data.M2);
set(ax,'YDir','normal','SortMethod','childorder'); hold(ax,'on'); colormap(ax,summer(256));
cb=colorbar(ax); cb.Label.String='M2'; r=front.ridge;
lambdaMax=max(abs(front.lambda_AT_minus)); ratio=min(1,abs(front.lambda_AT_minus(:))./lambdaMax);
gray=.75.*(1-ratio)+.03.*ratio;
for q=1:numel(r.x)-1
    plot(ax,r.x(q:q+1),r.y(q:q+1),'-','Color',repmat(mean(gray(q:q+1)),1,3), ...
        'LineWidth',2.8,'HandleVisibility','off');
end
colors={[0.68 0.68 0.68],[0.42 0.42 0.42],[0 0.25 0.95],[0.95 0.02 0.02]};
widths=[1.8 1.8 .9 1.2];
for j=1:numel(front.eigen_samples)
    s=front.eigen_samples(j); if any(~isfinite(s.AT_values)), continue; end
    plot(ax,s.x,s.y,'ko','MarkerFaceColor','w','MarkerSize',3,'HandleVisibility','off');
    vectors={s.S_vectors(:,1),s.S_vectors(:,2),s.AT_vectors(:,1),s.AT_vectors(:,2)};
    values=[s.S_values(1),s.S_values(2),s.AT_values(1),s.AT_values(2)];
    for q=1:4
        halfLength=maxHalfLengthKm*(abs(values(q))/globalMax)^(1/rootDegree);
        uv=halfLength.*vectors{q}; draw_arrow(ax,s.x,s.y,uv(1),uv(2),colors{q},widths(q),values(q)<0);
    end
end
margin=maxHalfLengthKm+3; axis(ax,'equal');
xlim(ax,[max(min(data.lonir(:)),min(r.x)-margin) min(max(data.lonir(:)),max(r.x)+margin)]);
ylim(ax,[max(min(data.latir(:)),min(r.y)-margin) min(max(data.latir(:)),max(r.y)+margin)]);
set(ax,'FontSize',11,'Box','on','Layer','top'); grid(ax,'on'); xlabel(ax,'x (km)'); ylabel(ax,'y (km)');
title(ax,sprintf('Retained FL: FL%d, %.2f km',front.fl_id,r.length_km));
exportgraphics(fig,outputFile,'Resolution',600,'BackgroundColor','white');
end

function draw_arrow(ax,x,y,ux,uy,c,w,inward)
if inward
    quiver(ax,x+ux,y+uy,-ux,-uy,0,'Color',c,'LineWidth',w,'MaxHeadSize',.5,'AutoScale','off','HandleVisibility','off');
    quiver(ax,x-ux,y-uy,ux,uy,0,'Color',c,'LineWidth',w,'MaxHeadSize',.5,'AutoScale','off','HandleVisibility','off');
else
    quiver(ax,x,y,ux,uy,0,'Color',c,'LineWidth',w,'MaxHeadSize',.5,'AutoScale','off','HandleVisibility','off');
    quiver(ax,x,y,-ux,-uy,0,'Color',c,'LineWidth',w,'MaxHeadSize',.5,'AutoScale','off','HandleVisibility','off');
end
end

