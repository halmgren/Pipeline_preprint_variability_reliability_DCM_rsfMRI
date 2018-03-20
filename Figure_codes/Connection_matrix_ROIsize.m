function Connection_matrix_ROIsize(SPM_dir,Work_dir)

for Size_ROI=[1 4]
    load([Work_dir '/Results_paper_variability/DCM/ROI_Size/Smith/Full_model/' num2str(Size_ROI*4) '/PEB_group/PEB_A_mean_comp_ROIsize_group_DMN.mat'],'PEB');
    load([Work_dir '/DatasetKuehn/sub-01_summary/DCM/Basic/Smith/Full_model/GCM_DMN_full_estim.mat']);
    DCM=GCM;
    
    s = get(0,'screensize');
    
    w = round(s(3) * 0.25);
    h = round(s(4) * 0.38);
    
    p = get(gcf,'OuterPosition');
    b = p(4) - h;
    
    l = p(1) - w;                 % Attempt left of figure window
    if l < 1, l = p(1)+p(3); end  % Attempt right of figure window
    if l > s(3), l = p(1); end    % Align with figure window
    
    xPEB.mtx_fig = figure('Name','Connectivity','NumberTitle','off','units','normalized',...
        'OuterPosition',[0 0 1 1],'Color','white',...
        'ToolBar','none');
    
    ci=spm_invNcdf(1-0.05);
    EP=full(vec2mat(PEB.Ep(1:16),4)');
    CP=diag(PEB.Cp);
    CP=full(vec2mat(CP(1:16),4)');
    sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));
    
    A_matrix=full(vec2mat(PEB.Ep(1:16),4)');
    
    A_matrix(sgn==-1)=NaN;
    
    h=imagesc(A_matrix);
    set(h,'alphadata',~isnan(A_matrix))
    
    colorbar;
    
    axis square;
    xlabel('\textbf{\underline{From}}','FontSize',40,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',40,'Fontweight','bold','Interpreter','latex');
set(gca,'XAxisLocation','top','Tag','connectivity');

    title_str = ['Connectivity'];
    
    title(title_str,'FontSize',46);
    
    
    for region=1:length({DCM{1}.xY.name})
        regions(region)={DCM{1}.xY(region).name(5:end)};
    end
    
    if size(A_matrix,1) == size(A_matrix,2) && ...
            size(A_matrix,1) == length({DCM{1}.xY.name})
        set(gca,'YTickLabel',regions,...
            'YTick',1:length({DCM{1}.xY.name}),...
            'XTickLabel',regions,'fontweight','bold','fontsize',40,'XTick',...
            1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
    end
    
    for side1=1:4
        for side2=1:4
            if ~isnan(A_matrix(side2,side1))
                text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
                
            elseif isnan(A_matrix(side2,side1))
                text(side1,side2,'ns','HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
            end
        end
    end
    cmap_custom=parula;
    cmap_custom(1:4,:)=[];
    colormap(cmap_custom);
    
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,[Work_dir '/Figures_paper_variability/group_ROI_size_' num2str(4*Size_ROI) '.bmp']);
    close;
    clear DCM PEB;
    
end

end