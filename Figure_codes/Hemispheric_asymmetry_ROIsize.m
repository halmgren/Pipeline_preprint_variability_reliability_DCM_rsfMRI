function Hemispheric_asymmetry_ROIsize(SPM_dir,Work_dir)

%%%%%%%%%%%%%
%group-lateral
%%%%%%%%%%%%%
for Size_ROI=[1 4]
    load([Work_dir '/Results_paper_variability/DCM/ROI_Size/Smith/Full_model/' num2str(Size_ROI*4)  '/PEB_group/Lateralization_index_A_comp_ROIsize_DMN.mat']);
    c=mean_diff;
    v=var_of_sum;
    PP=posterior_probability;
    
    T=0;
    x    = c + [-512:512]*sqrt(v)*6/512;
    p    = full(1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v)));
    i    = 1:1025;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(x,p,[1 1]*T,[0 max(p)],'-.k','linewidth',5);
    
    if c>0
        str  = sprintf('P(contrast > %0.2f) = %.1f%s',T,PP*100,'%');
    elseif c<0
        str  = sprintf('P(contrast < %0.2f) = %.1f%s',T,(1-PP)*100,'%');
    end
    
    title(str,'FontSize',24);
    xlabel('Contrast of parameter estimates');
    ylabel('Posterior density');
    xlim([-0.1 0.4]);
    hold on;
    fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*0.5);
    plot(x,p,[1 1]*T,[0 max(p)],'-.k','linewidth',5','color','k');
    axis square, grid on
    hold off;
    set(gca,'fontweight','bold','fontsize',36);
    
    %axes
    ax=gca;
    yrule=ax.YAxis;
    yrule.FontSize=44;
    xrule=ax.XAxis;
    xrule.FontSize=36;

    ylim([0 16]);
    
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,[Work_dir '/Figures_paper_variability/group_ROI_size_' num2str(4*Size_ROI) '_lateralization.bmp']);
    close;
end

end