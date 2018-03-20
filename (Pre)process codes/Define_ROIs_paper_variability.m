function [ROI_list]=Define_ROIs_paper_variability(name_ROI_def)

%naming is very important: first three letters should correspond to
%network, further letters correspond to region
%regions should contain three or four letter, no more, no less
if strcmp(name_ROI_def,'Smith')
    
    ROI_list={'DMN_PRC',[2,-58,30];...      %%Precuneus/posteriorPCC
        'DMN_mPFC',[2,56,-4];...
        'DMN_lIPC',[-44,-60,24];...    %inferior parietal
        'DMN_rIPC',[54,-62,28]};        %inferuir parietal
end

end