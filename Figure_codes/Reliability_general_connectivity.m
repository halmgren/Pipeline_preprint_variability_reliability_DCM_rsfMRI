function Reliability_general_connectivity(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%
%All permutations: Normal vs delateralized
%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/WS_correlation_DMN.mat']);
        cr_group(1,tmp)=mean(cr_6);
        
    end
end
cr_group(:,end+1)=mean(cr_group,2);

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/WS_correlation_delat_DMN.mat']);
        cr_group_delat(1,tmp)=mean(cr_6);
    end
end

cr_group_delat(:,end+1)=mean(cr_group_delat,2);

cr_group_both=[cr_group;cr_group_delat];

figure('units','normalized','outerposition',[0 0 1 1]);
a=plot(([cr_group(:,1:17);cr_group_delat(:,1:17)]),'-x','LineWidth',4,'MarkerSize',30);
hold on;
b=plot(([cr_group(:,18);cr_group_delat(:,18)]),'-xk','LineWidth',10,'MarkerSize',60);

leg_loc=legend([b],{'AVERAGE'});
set(leg_loc,'Position',get(leg_loc,'Position')-[0.15 +0.05 0 0])
set(gca,'XTick',[1 2],'YLim',[0 1],'XLim',[0.75 2.25],'fontsize',40,'fontweight','bold','XTickLabel',{'Original','Dominance-ordered'});
ylabel('Mean correlation','fontweight','bold','fontsize',44);

set(gcf,'PaperPositionMode','auto');
grid on;
saveas(gcf,[Work_dir '/Figures_paper_variability/Reliability_general_connectivity.bmp']);
close;
clear cr_group_delat cr_group_delat cr_group_both;

end