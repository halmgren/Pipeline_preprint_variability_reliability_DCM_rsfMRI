function Comparison_GSR_estimates(SPM_dir,Work_dir)

clear PEB;
procedure='Basic';
load([Work_dir '/Results_paper_variability/DCM/' procedure '/Smith/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_DMN.mat'],'PEB');
PEB1=PEB;
clear PEB;
procedure='GSR';
load([Work_dir '/Results_paper_variability/DCM/' procedure '/Smith/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_DMN.mat'],'PEB');
PEB2=PEB;

figure('units','normalized','outerposition',[0 0 1 1]);
spm_plot_ci_adapt_v3(PEB1.Ep,PEB1.Cp);
hold on;
spm_plot_ci_adapt_v2(PEB2.Ep,PEB2.Cp);
alpha 0.4;

xlabel('Connections','fontweight','bold','fontsize',55);
ylabel('Posterior expectations (Hz)','fontweight','bold','fontsize',32);

%axes
ax=gca;
yrule=ax.YAxis;
yrule.FontSize=32;
xrule=ax.XAxis;
xrule.FontSize=38;

set(gca,'Xtick',[1 5 10 15])
set(gcf,'PaperPositionMode','auto');

saveas(gcf,[Work_dir '/Figures_paper_variability/comparison_GSR.bmp']);
close;
end