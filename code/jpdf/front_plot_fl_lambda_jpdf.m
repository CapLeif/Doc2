function front_plot_fl_lambda_jpdf(stats,outputFile,kind,statistic)
%FRONT_PLOT_FL_LAMBDA_JPDF Plot FL distribution or unsmoothed angle conditionals.
fig=figure('Units','centimeters','Position',[2 2 18 16],'Color','w');
ax=axes(fig); c=stats.bin_centers;
switch lower(kind)
    case 'distribution'
        z=log10(stats.pdf); label='$\log_{10}\ PDF$';
        hm=imagesc(ax,c,c,z.');
        colormap(ax,gray(256));
        titleText=['Day 40 FL points on \lambda-JPDF (N = ' num2str(stats.sample_count) ')'];
    case 'theta'
        z=stats.(['theta_' lower(statistic)]); label=['$\theta$ conditional ' lower(statistic) ' ($^\circ$)'];
        hm=imagesc(ax,c,c,z.'); colormap(ax,nclCM(353,16)); caxis(ax,[0 180]);
        titleText=['Day 40 FL: conditional ' lower(statistic) ' of \theta'];
    case 'epsilon'
        z=stats.(['epsilon_' lower(statistic)]); label=['$\epsilon(V_-,\nabla_h b)$ conditional ' lower(statistic) ' ($^\circ$)'];
        hm=imagesc(ax,c,c,z.'); colormap(ax,nclCM(353,16)); caxis(ax,[0 90]);
        titleText=['Day 40 FL: conditional ' lower(statistic) ' of \epsilon'];
    otherwise
        error('Front:UnknownPlotKind','Unknown plot kind: %s',kind);
end
set(hm,'AlphaData',isfinite(z.')); set(ax,'YDir','normal','Color',[.92 .92 .92], ...
    'FontName','Arial','FontSize',12,'LineWidth',1.1,'Box','on'); hold(ax,'on');
lims=[stats.bin_edges(1) stats.bin_edges(end)];
plot(ax,lims,lims,'k--','LineWidth',0.9); plot(ax,lims,-lims,'k--','LineWidth',0.9);
xline(ax,0,'k--','LineWidth',0.9); yline(ax,0,'k--','LineWidth',0.9);
xlim(ax,lims); ylim(ax,lims); axis(ax,'square');
xlabel(ax,'$\Delta/f$','Interpreter','latex');
ylabel(ax,'$-\sqrt{\Omega}/|f|$','Interpreter','latex');
title(ax,titleText,'FontWeight','normal'); cb=colorbar(ax); cb.Label.Interpreter='latex'; cb.Label.String=label;
text(ax,-8.7,-4.0,'Cs-','FontWeight','bold'); text(ax,-2.3,-4.0,'Sc-','FontWeight','bold');
text(ax,1.1,-4.0,'Sd-','FontWeight','bold'); text(ax,7.5,-4.0,'Ds+','FontWeight','bold');
text(ax,.02,.02,sprintf('bin = %.2f; no smoothing',diff(stats.bin_edges(1:2))), ...
    'Units','normalized','FontSize',9,'BackgroundColor','w','Margin',2);
exportgraphics(fig,outputFile,'Resolution',600,'BackgroundColor','white'); close(fig);
end

