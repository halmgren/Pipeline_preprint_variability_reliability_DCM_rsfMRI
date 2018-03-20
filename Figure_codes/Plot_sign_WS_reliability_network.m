function Plot_sign_WS_reliability_network(dataset,subject,ntwrk_size,ntwrk_VOI_names,ntwrk_VOI_coord,network_name,name_ROI_def,stability_threshold,procedure,SPM_dir,Work_dir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adapted from spm_dcm_display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/Sign_stability_' network_name '.mat'],'Positive_sign_prop','Negative_sign_prop','Ns_sign_prop','Number_HQ_sessions');
          
a(:,:,1)=Positive_sign_prop;
a(:,:,2)=NaN(4);

clear DCM xY;
load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/GCM_' network_name '_full_estim.mat']);

for VOI1=1:ntwrk_size
    for VOI2=1:ntwrk_size
        if strcmp(ntwrk_VOI_names{VOI1},GCM{1}.xY(VOI2).name(5:end))
            DCM.xY(VOI2).name=ntwrk_VOI_names{VOI1};
            DCM.xY(VOI2).xyz=ntwrk_VOI_coord{VOI1}';
        end
    end
end

%If image is unclear, some coordinates are changed
DCM.xY(1).xyz=[2 -35 30]';

xY=DCM.xY;

n = 2;

ha = gca;

% graphics parameters
%--------------------------------------------------------------------------
col   = {'r','g','b','c','y','m'};
rad   = 20;                   % radius of self-connections
w     = 2;                    % line width
M     = 10;                   % MarkerSize for regions
Upper     = stability_threshold;                  % Display threshold on p-value

% get dimensions, locations and names
%--------------------------------------------------------------------------
m     = size(xY,2);
L     = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name(1:min(end,3));
end
L     = [L; ones(1,m)];
o     = mean(L,2);
M1    = spm_matrix(-o(1:3)');


% get orientation with the greatest dispersion
%--------------------------------------------------------------------------
for i = 1:3
    u{i}      = eye(4,3);
    u{i}(:,i) = [];
    s(i)      = det(u{i}'*M1*L*L'*M'*u{i});
end
[i,j]   = max(s);
u       = u{j};

% compute projection matrix for 'principal' plane
%--------------------------------------------------------------------------
M2      = u';
M2(4,4) = 1;
M1(j,4) = 0;
L       = M2*M1*L;

% coordinates
%--------------------------------------------------------------------------
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


% get T1 background
%--------------------------------------------------------------------------
V       = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
ijk     = V.mat \ xyz;
t1      = spm_sample_vol(V,ijk(1,:),ijk(2,:),ijk(3,:),2);
t1      = (64 - 16) + 16*t1/max(t1(:));

t1(find(t1==48))=64;

% Watermark and regions
%--------------------------------------------------------------------------
str     = get(get(ha,'Title'),'String');
image(rot90(reshape(t1,length(i),length(j))),'parent',ha);
colormap(gray);
axis(ha,'image','off')
title(ha,str)

% Connections
%--------------------------------------------------------------------------
Q     = [-pi:pi/32:pi];
Q     = rad*[sin(Q); cos(Q)];
q     = 1/3;
disp(a);
for i = 1:length(a)
    for j = 1:length(a)
        if ~isnan(a(i,j,1))
            
            % show connections
            %--------------------------------------------------------------
            if i ~= j
                
                % line
                %----------------------------------------------------------
                k = rem(j - 1,length(col)) + 1;
                
                
                % if significant
                %----------------------------------------------------------
                if (a(i,j,1) >= Upper)||(Negative_sign_prop(i,j) >= Upper)
                    
                    if i>j
                        h = line(L(1,[i j]),L(2,[i j]),'Color',col{k},...
                            'LineStyle',':',...
                            'LineWidth',w);
                        set(h,'LineStyle','-','LineWidth',w)
                    else
                        h = line(L(1,[i j])+2,L(2,[i j])+2,'Color',col{k},...
                            'LineStyle',':',...
                            'LineWidth',w);
                        set(h,'LineStyle','-','LineWidth',w)
                    end
                    
                    % text
                    %------------------------------------------------------
                    u     = 0.3*(L(1,j) - L(1,i)) + L(1,i);
                    v     = 0.3*(L(2,j) - L(2,i)) + L(2,i);
                    str   = {};
                    for k = 1
                        if a(i,j,1)>=Upper
                            str{k} = sprintf('+');
                            h     = text(u,v,1,str(:),'FontSize',16,'Fontweight','bold',...
                                'HorizontalAlignment','Center');
                        else
                            str{k} = sprintf('-');
                            h     = text(u,v,1,str(:),'FontSize',24,'Fontweight','bold',...
                                'HorizontalAlignment','Center');
                        end
                    end
                    
                end
                
            else
                
                % line
                %----------------------------------------------------------
                k     = rem(i - 1,length(col)) + 1;
                u     = L(1,i);
                v     = L(2,i);
                u     = Q(1,:) + u;
                v     = Q(2,:) + v;
                
                
            end
        end
    end
end


% projected coordinates of voxels within region[s]
%--------------------------------------------------------------------------
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
if exist('sbplt')
    title(sbplt,[dataset(8:end) '\_subj' num2str(subject)]);
end
if exist('sbplt2')
    title(sbplt2,[dataset(8:end) '\_subj' num2str(subject)]);
end
set(gcf,'PaperPositionMode','auto');

saveas(gcf,[Work_dir '/Figures_paper_variability/Sign_reliability_network_' dataset '_' num2str(subject) '.bmp']);


end