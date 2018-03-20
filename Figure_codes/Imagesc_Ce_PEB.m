function Imagesc_Ce_PEB(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Group between subject variability (PEB.Ce)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/PEB_A_mean_group_DMN.mat'],'PEB');
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

xPEB.mtx_fig = figure('Name','Connectivity','NumberTitle','off',...
    'OuterPosition',[l b w h],'Color','white',...
    'ToolBar','none');

% ci=spm_invNcdf(1-0.05);
% EP=full(vec2mat(PEB.Ep(1:16),4)');
% CP=diag(PEB.Cp);
% CP=full(vec2mat(CP(1:16),4)');
% sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));

Ce_matrix=full(vec2mat(diag(PEB.Ce),4)');

% Ce_matrix(sgn==-1)=NaN;
figure('units','normalized','outerposition',[0 0 1 1]);
h=imagesc(Ce_matrix*100,[0.2 6.9]);
set(h,'alphadata',~isnan(Ce_matrix))

colorbar;

axis square;
xlabel('\textbf{\underline{From}}','FontSize',40,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',40,'Fontweight','bold','Interpreter','latex');
set(gca,'XAxisLocation','top','Tag','connectivity');

title_str = ['Between-subject variance (x10^{-2})'];

title(title_str,'FontSize',46);

%         if size(A_matrix,1) == size(A_matrix,2) && ...
%                 size(A_matrix,1) == length({DCM{1}.xY.name})
%             set(gca,'YTickLabel',{DCM{1}.xY.name},...
%                 'YTick',1:length({DCM{1}.xY.name}),...
%                 'XTickLabel',{''},'TickLabelInterpreter', 'none');
%         end

for region=1:length({DCM{1}.xY.name})
    regions(region)={DCM{1}.xY(region).name(5:end)};
end

if size(Ce_matrix,1) == size(Ce_matrix,2) && ...
        size(Ce_matrix,1) == length({DCM{1}.xY.name})
    set(gca,'YTickLabel',regions,...
        'YTick',1:length({DCM{1}.xY.name}),...
        'XTickLabel',regions,'fontweight','bold','fontsize',40,'XTick',...
        1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
end

for side1=1:4
    for side2=1:4
        if ~isnan(Ce_matrix(side2,side1))
            text(side1,side2,num2str(round(Ce_matrix(side2,side1)*100,2)),'HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
            
        elseif isnan(Ce_matrix(side2,side1))
            text(side1,side2,'ns','HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
        end
    end
end
cmap_custom=parula;
cmap_custom(1:4,:)=[];
colormap(cmap_custom);

set(gcf,'PaperPositionMode','auto');
saveas(gcf,[Work_dir '/Figures_paper_variability/Group_PEB_Ce.bmp']);
close;
clear DCM PEB;

end