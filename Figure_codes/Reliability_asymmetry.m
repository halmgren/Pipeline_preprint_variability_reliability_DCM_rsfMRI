function Reliability_asymmetry(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function uses gramm-toolbox
%Downloaded from https://github.com/piermorel/gramm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lateralization ordered based on stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=0;
tmp2=-1;

for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        tmp2=tmp2+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These specific subjects were excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            tmp2=tmp2-1;
            continue
        end
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/Lateral_index_individ_DMN.mat']); %load lateralization
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']); %load QC's
        
        
        flag=zeros(1,length(mean_diff));
        for diagn=1:length(mean_diff)
            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                flag(diagn)=1;
            end
            
        end
        
        mean_diff(flag==1)=[];
        posterior_probability(flag==1)=[];
        flag2=zeros(1,length(mean_diff));
        
        for diagn=1:length(mean_diff)
            if posterior_probability(diagn)<0.95&&posterior_probability(diagn)>0.05
                flag2(diagn)=1;
            end
        end
        
        mean_diff_group{tmp}=mean_diff;
        stab(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff);
        Non_sign{tmp,:}=flag2;
    end
end

[~,order]=sort(stab);
for i=1:length(stab)
    Non_sign_order{i}=Non_sign{order(i)};
    mean_diff_group_order{i}=mean_diff_group{order(i)};
    
    stab2(i)=stab(order(i));
end

clear x_axis y_axis A g;

tmp=0;
for subject=1:17
    tmp=tmp+1;
    
    x_axis(tmp:tmp+length(mean_diff_group_order{subject})-1)=ones(size(mean_diff_group_order{subject}))*subject;
    tmp=tmp+length(mean_diff_group_order{subject})-1;
end
x_axis=x_axis';

tmp=0;
for subject=1:17
    
    tmp=tmp+1;
    if rem(subject,2)
        clr=1;
    else
        clr=2;
    end
    clrs(tmp:tmp+length(mean_diff_group_order{subject})-1)=ones(size(mean_diff_group_order{subject}))*clr;
    tmp=tmp+length(mean_diff_group_order{subject})-1;
end
clrs=clrs';

tmp2=0;
for subject=1:17
    
    for session=1:length(Non_sign_order{subject})
        tmp2=tmp2+1;
        if Non_sign_order{subject}(session)==1
            clrs(tmp2)=3;

        end
    end
end

tmp=0;
for subject=1:17
    tmp=tmp+1;
    
    y_axis(tmp:tmp+length(mean_diff_group_order{subject})-1)=mean_diff_group_order{subject};
    tmp=tmp+length(mean_diff_group_order{subject})-1;
end
y_axis=y_axis';

g(1,1)=gramm('x',x_axis,'y',y_axis,'color',clrs,'size',ones(1,length(x_axis))*2);
fig1=figure('units','normalized','outerposition',[0 0 1 1]);
g(1,1).geom_hline('yintercept',0,'style','-k');
g(1,1).geom_jitter('width',0.6,'dodge',0.5);
g(1,1).set_names('x','Subject (ordered)','y','Asymmetry index');
g(1,1).set_color_options('map',[0 0 0;0.75 0.75 0.75;0 0.5 1]);
g(1,1).no_legend;
g(1,1).set_point_options('base_size',10);
g(1,1).axe_property('FontWeight','bold','FontSize',36)
g(1,1).set_layout_options('Position',[0.1 0.1 0.8 0.8],'margin_height',[0.15 0.02],'margin_width',[0.1 0.02],'legend_pos',[0 0 0 0],'redraw',false);
a=g(1,1).draw;
A=a.facet_axes_handles;


Ylim=A.YLim;
set(A,'YLim',[Ylim(1), Ylim(2)+0.7])

for i=1:17
    line('XData',[i i],'Ydata',[Ylim(2)+0.2 Ylim(2)+0.2+length(mean_diff_group_order{i})/(157/0.5)],'LineWidth',12,'color','k','parent',A);
end
set(A,'Ytick',[-1.5:0.5:1.5]);

set(gcf,'PaperPositionMode','auto');
saveas(gcf,[Work_dir '/Figures_paper_variability/WS_lateralization_reliability.bmp']);
close;

end