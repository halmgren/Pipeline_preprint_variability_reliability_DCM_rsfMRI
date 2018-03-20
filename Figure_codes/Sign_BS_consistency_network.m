function Sign_BS_consistency_network(prop_posi,prop_nega,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Figure for sign_stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a(:,:,1)=prop_posi;
a(:,:,2)=prop_nega;

load([Work_dir '/DatasetKuehn/sub-01/ses-001/func/DCM/Basic/Smith/Full_model/GCM_full_estim.mat']); %Used to determine location of regions


[ROI_list]=Define_ROIs_paper_variability('Smith');

tmp=0;
for VOI_number=1:size(ROI_list,1)
    ntwrk=ROI_list{VOI_number,1}(1:3);
    
    if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
        n=n+1;
        ntwrk_size(tmp)=ntwrk_size(tmp)+1;
        ntwrk_VOI_names{n,tmp}=ROI_list{VOI_number,1}(5:end);
        ntwrk_VOI_coord{n,tmp}=ROI_list{VOI_number,2};
        continue
        
    else
        n=1;
        tmp=tmp+1;
        ntwrk_size(tmp)=1;
        ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
        ntwrk_VOI_names{n,tmp}=ROI_list{VOI_number,1}(5:end);
        ntwrk_VOI_coord{n,tmp}=ROI_list{VOI_number,2};
    end
end

network_number=1;
for VOI1=1:ntwrk_size(network_number)
    for VOI2=1:ntwrk_size(network_number)
        if strcmp(ntwrk_VOI_names{VOI1,network_number},GCM{1}.xY(VOI2).name(5:end))
            DCM.xY(VOI2).name=ntwrk_VOI_names{VOI1,network_number};
            DCM.xY(VOI2).xyz=ntwrk_VOI_coord{VOI1,network_number}';
        end
    end
end

DCM.xY(1).xyz=[2 -35 30]'; %Make picture more clear

n = 2;

figure;


ha = gca;

% graphics parameters
%--------------------------------------------------------------------------
col   = {'r','g','b','c','y','m'};
rad   = 20;                   % radius of self-connections
w     = 2;                    % line width
M     = 10;    

xY=DCM.xY;

xY(2).name='dIP';
xY(3).name='mPF';
xY(4).name='ndIP';
m     = size(xY,2);
L     = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name(1:min(end,5));
end
L     = [L; ones(1,m)];
o     = mean(L,2);
M1    = spm_matrix(-o(1:3)');

for i = 1:3
    u{i}      = eye(4,3);
    u{i}(:,i) = [];
    s(i)      = det(u{i}'*M1*L*L'*M'*u{i});
end
[i,j]   = max(s);
u       = u{j};

M2      = u';
M2(4,4) = 1;
M1(j,4) = 0;
L       = M2*M1*L;

i       = (min(L(1,:)) - M):(max(L(1,:)) + M);
j       = (min(L(2,:)) - M):(max(L(2,:)) + M);
i=-90:90;
j=-90:125;
x       = kron(ones(size(j)),i);
y       = kron(j,ones(size(i)));
xyz     = [x; y; zeros(1,length(x)); ones(1,length(x))];
xyz     = pinv(M1)*pinv(M2)*xyz;
M3      = [1 0 0 -min(i); 0 -1 0 max(j); 0 0 0 0; 0 0 0 1];
L       = M3*L;

V       = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
ijk     = V.mat \ xyz;
t1      = spm_sample_vol(V,ijk(1,:),ijk(2,:),ijk(3,:),2);
t1      = (64 - 16) + 16*t1/max(t1(:));


t1(find(t1==48))=64;

% Watermark and regions
%--------------------------------------------------------------------------
str     = get(get(ha,'Title'),'String');
image(rot90(reshape(t1,length(i),length(j))),'parent',ha)%
colormap(gray);
axis(gca,'image','off');
title(ha,str)

Q     = [-pi:pi/32:pi];
Q     = rad*[sin(Q); cos(Q)];
q     = 1/3;

for i = 1:length(a)
    for j = 1:length(a)
        if ~isnan(a(i,j,1))
%             pause;
            hold on;
            % show connections
            %--------------------------------------------------------------
            if i ~= j
                
                % line
                %----------------------------------------------------------
                k = rem(j - 1,length(col)) + 1;
                
                
                % if significant
                %----------------------------------------------------------
                if (a(i,j,1) >= 0.60)||(a(i,j,2) >= 0.60)
                    
                    %                     if i>j
                    %                         if a(i,j,1)>0
                    %                             clr='g';
                    %                         else
                    %                             clr='r';
                    %                         end
                    
                    if a(i,j,1) >= 0.75
                        colo='g';
                        
                        h = line(L(1,[i j]),L(2,[i j]),'Color',col{k},...
                            'LineStyle',':',...
                            'LineWidth',w);
                        set(h,'LineStyle','-','LineWidth',w)
                        
                     elseif a(i,j,2) >= 0.75
                        colo='r';
                        
                        h = line(L(1,[i j]+1),L(2,[i j]),'Color',col{k},...
                            'LineStyle',':',...
                            'LineWidth',w);
                        set(h,'LineStyle','-','LineWidth',w)
                    elseif (a(i,j,1) >= 0.60)
                        colo='g';
                        h = line(L(1,[i j])+2,L(2,[i j])+2,'Color',col{k},...
                            'LineStyle','--',...
                            'LineWidth',w);
                        set(h,'LineStyle','--','LineWidth',w)
                        
                    elseif a(i,j,2) >= 0.60
                        colo='r';
                           h = line(L(1,[i j])+2,L(2,[i j])+2,'Color',col{k},...
                            'LineStyle','--',...
                            'LineWidth',w);
                        set(h,'LineStyle','--','LineWidth',w)
                         
                    end
                end
                
                
                % text
                %------------------------------------------------------
                u     = 0.17*(L(1,j) - L(1,i)) + L(1,i);
                v     = 0.17*(L(2,j) - L(2,i)) + L(2,i);
                str   = {};
                for k = 1%:size(a,3)
                    if a(i,j,1)>=0.60
                        str{k} = sprintf('+');
                        
                        if (i==3)&&(j==2)
                            
                            h     = text(u-4,v,1,str(:),'FontSize',20,'Fontweight','bold',...
                            'HorizontalAlignment','Center');
                        elseif (i==3)&&(j==1)
                           h     = text(u+5,v,1,str(:),'FontSize',20,'Fontweight','bold',...
                            'HorizontalAlignment','Center');
                        elseif (i==1)&&(j==2)
                           h     = text(u-8,v+1,1,str(:),'FontSize',20,'Fontweight','bold',...
                            'HorizontalAlignment','Center');
                        elseif (i==4)&&(j==2)
                           h     = text(u-4,v-4,1,str(:),'FontSize',20,'Fontweight','bold',...
                            'HorizontalAlignment','Center');
                        else
                            
                            h     = text(u,v,1,str(:),'FontSize',20,'Fontweight','bold',...
                            'HorizontalAlignment','Center');
                        end
                    elseif a(i,j,2)>=0.60
                        str{k} = sprintf('-');
                        
                        h     = text(u+2,v+4,1,str(:),'FontSize',28,'Fontweight','bold',...
                            'HorizontalAlignment','left');
                    end
                end
                %                 end
                
            end
        end
    end
end
set(gca,'Ydir','reverse');

hold on
for i = 1:m
    k = rem(i - 1,length(col)) + 1;
    line(L(1,i),L(2,i),...
        'Color',col{k},...
        'Marker','.',...
        'LineStyle','none',...
        'MarkerSize',98);
    
    text(L(1,i),L(2,i),name{i},'FontSize',16,...
        'FontWeight','Bold',...
        'Color','k',...
        'HorizontalAlignment','center')
end
hold off

set(gcf,'PaperPositionMode','auto');
saveas(gcf,[Work_dir '/Figures_paper_variability/BS_sgn_consistency_network.bmp']);
close;

end

