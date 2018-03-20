function [prop_posi prop_nega]=Sign_BS_consistency(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count positives and negatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure='Basic';
name_ROI_def='Smith';
ntwrk_name='DMN';

clear A_matrix;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name '.mat']);
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name '.mat']);
        catch ME
            tmp=tmp-1;
            disp('PEB group')
            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
            %
            continue;
            
        end
        %We don't include participants with less then 8
        %useful sessions.
        if length(PEB.Snames)<8
            tmp=tmp-1;
            continue
        else 
            
            PEB_group{tmp}=PEB;
            Lat_group{tmp}=mean_diff;
        end
        A_matrix_tmp=full(vec2mat(PEB.Ep(1:16),4)');
        
        %CI
        ci=spm_invNcdf(1-0.05);
        EP=full(vec2mat(PEB.Ep(1:16),4)');
        CP=diag(PEB.Cp);
        CP=full(vec2mat(CP(1:16),4)');
        sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));
        
        A_matrix_tmp(sgn==-1)=NaN;
        A_matrix(:,:,tmp)=A_matrix_tmp;
        clear ci EP CP sgn;
    end
end

A_matrix_delat=A_matrix;


for subject=1:length(Lat_group)
    if Lat_group{subject}>0
        A_matrix_delat([2 4],:,subject)=A_matrix_delat([4 2],:,subject);
        A_matrix_delat(:,[2 4],subject)=A_matrix_delat(:,[4 2],subject);
    end
end

for sz1=1:4
    for sz2=1:4
        number_nega(sz1,sz2)=sum(squeeze(sign(A_matrix_delat(sz1,sz2,:)))==-1);
        number_posi(sz1,sz2)=sum(squeeze(sign(A_matrix_delat(sz1,sz2,:)))==1);
        number_NaN(sz1,sz2)=sum(isnan(squeeze(A_matrix_delat(sz1,sz2,:))));
        
    end
end

%Make second column and row dominant and fourth non-dominant
number_nega([2 4],:)=number_nega([4 2],:);
number_nega(:,[2 4])=number_nega(:,[4 2]);
number_posi([2 4],:)=number_posi([4 2],:);
number_posi(:,[2 4])=number_posi(:,[4 2]);
number_NaN([2 4],:)=number_NaN([4 2],:);
number_NaN(:,[2 4])=number_NaN(:,[4 2]);

%Ignore self-connections
number_posi(logical(eye(size(number_posi))))=0;
number_nega(logical(eye(size(number_nega))))=0;
number_NaN(logical(eye(size(number_NaN))))=0;

%proportion iso counts
prop_posi=number_posi/tmp;
prop_nega=number_nega/tmp;
prop_NaN=number_NaN/tmp;



cmap=[1 1 1; 0 1/2 0; 0 1 0; 1/2 0 0; 1 0 0];

sign_stab=zeros(4,4);
for sz1=1:4
    for sz2=1:4
        if prop_posi(sz1,sz2)>0.8
            sign_stab(sz1,sz2)=1;
        elseif prop_posi(sz1,sz2)>0.6
            sign_stab(sz1,sz2)=2;
        elseif prop_nega(sz1,sz2)>0.8
            sign_stab(sz1,sz2)=3;
        elseif prop_nega(sz1,sz2)>0.6
            sign_stab(sz1,sz2)=4;
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);
h=imagesc(sign_stab);
colormap(cmap);

axis square;
xlabel('\textbf{\underline{From}}','FontSize',40,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',36,'Fontweight','bold','Interpreter','latex');
set(gca,'XAxisLocation','top','Tag','connectivity');

title_str = ['Connectivity'];

title(title_str,'FontSize',46);

for region=1:length({DCM{1}.xY.name})
    regions(region)={DCM{1}.xY(region).name(5:end)};
end

regions{2}='dIPC';
regions{4}='ndIPC';
if size(A_matrix,1) == size(A_matrix,2) && ...
        size(A_matrix,1) == length({DCM{1}.xY.name})
    set(gca,'YTickLabel',regions,...
        'YTick',1:length({DCM{1}.xY.name}),...
        'XTickLabel',regions,'fontweight','bold','fontsize',36,'XTick',...
        1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
end

for side1=1:4
    for side2=1:4
        if side1==side2
            continue
        end
        if ~isnan(A_matrix(side2,side1))
            text(side1,side2,[num2str(number_posi(side2,side1)) '(+)'],'HorizontalAlignment','center','Fontweight','bold','Fontsize',36,'VerticalAlignment','bottom');
        end
        
        if ~isnan(A_matrix(side2,side1))
            text(side1,side2,[num2str(number_nega(side2,side1)) '(-)'],'HorizontalAlignment','center','Fontweight','bold','Fontsize',36,'VerticalAlignment','top');
            
        end
        
    end
end

set(gcf,'PaperPositionMode','auto');
saveas(gcf,[Work_dir '/Figures_paper_variability/BS_sign_consistency.bmp']);
close;
clear DCM PEB;

end