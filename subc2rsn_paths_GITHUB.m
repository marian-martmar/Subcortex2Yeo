function P = subc2rsn_paths()
% Single declaration of paths to any folder/external file of the pipeline subc2rsn_...
% -------------------------------------------------------------------------
% Every pipeline script starts with  P = subc2rsn_paths();  and only uses P.<field>.
% TODO: If it needs GENPATH delcaration, addpath(genpath('P.<field>'), indicate with P.****_gp
% The comment after each field says which script(s) use it.

% To run the pipeline on another machine, EDIT THIS FILE ONLY
% -------------------------------------------------------------------------

P = struct();

%% ---- Toolboxes / external functions ------------------------------------
P.repo_functions    = '...\...\0.Functions\1.subc2rsn';                         % General folder with functions of the project
P.spm12             = '...\...\0.Functions\spm12\';                                                    % step 0
P.functions_tian    = '...\...\0.Functions\1.subc2rsn\functions_tian\';         % steps 0, 1 (show_progress, icaDim in /wishart)
P.functions_wishart = '...\...\0.Functions\1.subc2rsn\functions_tian\wishart\'; % step 0

%% ---- Atlases ------------------------------------------------------------
P.atlas_rsn     = '...\...\0.Data\1.Atlases\RSN';            % step 0     (Schaefer2018 1000 parcels, 7 and 17 networks, + centroid CSVs)
P.atlas_dk84    = '...\...\0.Data\1.Atlases\DK84';           % steps 0,3,4,5 (IIT_GM_Desikan_atlas.nii.gz, LUT_GM_Desikan_0to1.txt)
P.atlas_dbs80   = '...\...\0.Data\1.Atlases\DBS80';          % steps 0,3,4,5 (dbs80symm_2mm.nii, dbs80_nodes.xlsx)
P.atlas_tian    = '...\...\0.Data\1.Atlases\Tian_Github';    % steps 0,3,4,5 (Tian_Subcortex_S{1..4}_3T.nii + _label.txt)
P.atlas_aal     = '...\...\0.Data\1.Atlases\aal_for_SPM12\'; % steps 0,3

%% ---- HCP data (100 unrelated, MNI volume, 2 mm) ------------------------
P.hcp_local    = '...\...\0.Data\1.HCP_100unrelated_volume\'; % step 0 (SubjectIDs_100.mat array with HCP 100 unrelated IDs)
P.hcp_external = '...\...\HCP_100unrelated_volume\';          % step 1 (<ID>_rfMRI_REST1_LR_hp2000_clean.nii.gz volumetric timeseries)

%% ---- Step 1 inputs / outputs (trailing separator is REQUIRED: step 1 builds file names by concatenation)
P.outputs_0 = '...\...\1.Subc2rsn\outputs_0\';   % Input for step 1
P.dn_out             = '...\...\Filtered_dn\';   % step 1 (folder is created in case smoothed-denoised are saved; NOT FOR NOW)
P.dn_f_out           = '...\...\Filtered_dn_f\'; % step 1 (<ID>_xRL_dn_f.mat folder where we save smoothed-denoised+time filtered)
P.outputs_1 = '...\...\1.Subc2rsn\outputs_1\';   % Outputes of step 1
%% ---- Step 2 (CLUSTER) --------------------------------------------------
P.cluster_functions = '.../subc2rsn/Functions/'
P.cluster_ts_dir    = '.../Data/HCP_100unrelated_volume_dn_f/'; % step 2 (load xRL_dn_f files + Variables_for_FC_...mat)
P.cluster_out       = '.../subc2rsn/Outputs/';                  % step 2 where we save outputs: individual/group-avg FC(subc x rsn), rsn label assinged(subc x 1) % OLD Tian_DBS80_AAL_DK84/FC_dn_f_Results/FC_dn_f_Matrices/'; 

%% ---- Steps 3-8 (local, after copying the cluster results back) ----------
P.results_cluster= '...\...\1.Subc2rsn\outputs_2_cluster\';      % step 3
P.fc_dir         = '...\...\where_FC_computed_in_cluster';       % steps 3, 5 (cumFC_*.mat, <ID>_avgFC7/17.mat copied from the cluster)
P.tables_dir     = '...\...\1.Subc2rsn\Results_Tables_Figures\'; % parent folder of all results of steps 3-8

% Output stage (steps 3-8): ONE folder for results, tables and figures
P.results_dir = P.tables_dir;                              % subc2rsn_results.mat, subc2rsn_bootstrap_CI.mat, subc2rsn_views.mat
P.out_tables  = fullfile(P.results_dir, 'tables');         % step 4 and 5 (.csv, Map_*.mat)
P.out_figures = fullfile(P.results_dir, 'figures');        % steps 6, 7, 8 (.svg)
P.keyfig_outDir   = '...\...\OLD\SVG_Tian_KeyFigures';     % step 7a, MODE 'from_fig' (the .fig files of the old interactive script, to recover the cameras)
end
