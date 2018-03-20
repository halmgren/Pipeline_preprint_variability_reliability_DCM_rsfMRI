function Figures_paper_variability(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures paper 'Within and Between subject variability of the 'core' DMN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir([Work_dir '/Figures_paper_variability']);

%%%%%%%%%%%
%Figure 1
%%%%%%%%%%%
%Left panel
Mesh_DMN(SPM_dir,Work_dir);

%Middle panel
Imagesc_group_PEB(SPM_dir,Work_dir);

%Right panel;
Imagesc_Ce_PEB(SPM_dir,Work_dir);
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2: (continuation: see R-scripts)
%%%%%%%%%%%%%%%%%%%%%%%%%
Copy_for_R_script(SPM_dir,Work_dir);

%%%%%%%%%%%
%Figure 3
%%%%%%%%%%%
BS_violin_plots(SPM_dir,Work_dir);

%%%%%%%%%%%
%Figure 4
%%%%%%%%%%%
%left panel:
[prop_posi, prop_nega]=Sign_BS_consistency(SPM_dir,Work_dir);

%right panel:
Sign_BS_consistency_network(prop_posi,prop_nega,SPM_dir,Work_dir)

%%%%%%%%%%%%
%Figure 5
%%%%%%%%%%%%
%Upper panel:
Reliability_asymmetry(SPM_dir,Work_dir)

%Lower panel(s):
Reliability_asymmetry_extremes(SPM_dir,Work_dir)

%%%%%%%%%%
%Figure 6
%%%%%%%%%%
Reliability_general_connectivity(SPM_dir,Work_dir)

%%%%%%%%%%%%
%Figure 7
%%%%%%%%%%%%
Sign_WS_reliability_network(SPM_dir,Work_dir)

Plot_largest_connection(SPM_dir,Work_dir)

%%%%%%%%%%%%
%Figure 8
%%%%%%%%%%%%
Connection_matrix_ROIsize(SPM_dir,Work_dir)

Hemispheric_asymmetry_ROIsize(SPM_dir,Work_dir)
close all;
%%%%%%%%%%%
%Figure 9
%%%%%%%%%%%
Comparison_GSR_estimates(SPM_dir,Work_dir)


end