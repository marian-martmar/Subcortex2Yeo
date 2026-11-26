function C = subc2rsn_config()
% Constants shared by the output stage (steps 3 to 8 of pipeline subc2rsn_...)
% -------------------------------------------------------------------------
% Set only once: network names and colours, which group estimator the outputs show,
% file names of results, figure style...
%
%   C = subc2rsn_config();      % also adds <repository>/functions to the path
%   C.P contains machine-specific folders defined in subc2rsn_paths.m
% -------------------------------------------------------------------------

C = struct();
C.P = subc2rsn_paths();

%% ---- Networks: names (order = column order of every FC / count matrix) ----
%       Appeareance order in Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS
C.rsn7.labels  = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
%       Appeareance order in Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS
C.rsn17.labels = {'VisCent','VisPeri','SomMotA','SomMotB','DorsAttnA','DorsAttnB','SalVentAttnA','SalVentAttnB',...
                  'LimbicB','LimbicA','ContA','ContB','ContC','DefaultA','DefaultB','DefaultC','TempPar'};

%% ---- Networks: colours (rows in the SAME order as the labels) -------------
% Yeo-7: the table used in all the original figure scripts.
C.rsn7.colors = [ ...
 120  18 134;       % Vis
  70 130 180;       % SomMot
   0 118  14;       % DorsAttn
 196  58 250;       % SalVentAttn
 220 248 164;       % Limbic
 230 148  34;       % Cont
 205  62  78] / 255;% Default
% Yeo-17 in the SCHAEFER order of C.rsn17.labels using Schaefer2018 colour table
% (CBIG repository, MNI/freeview_lut), checked on 16/11/2026. Old script useed the original 17-network script listed these colours in the FreeSurfer Yeo2011
C.rsn17.colors = [ ...
 120  18 134;       % VisCent
 255   0   0;       % VisPeri
  70 130 180;       % SomMotA
  42 204 164;       % SomMotB
  74 155  60;       % DorsAttnA
   0 118  14;       % DorsAttnB
 196  58 250;       % SalVentAttnA
 255 152 213;       % SalVentAttnB
 122 135  50;       % LimbicB
 220 248 164;       % LimbicA
 230 148  34;       % ContA
 135  50  74;       % ContB
 119 140 176;       % ContC
 255 255   0;       % DefaultA
 205  62  78;       % DefaultB
   0   0 130;       % DefaultC
  12  48 255] / 255;   % TempPar

%% ---- Which group estimator the tables and figures show --------------------
% 'Z'  = Fisher average across subjects   (avgZ,  as7_Z,  as17_Z)   old tests used it for: exported maps, pies, S1 voxel figures
% 'FC' = arithmetic average of r          (avgFC, as7_FC, as17_FC)  old tests used it for: alluvial flows, S4 voxel figures and bootstrap CIs
C.estimator = 'Z'; % One switch for ALL outputs

%% ---- Subcortical Atlases ---------------------------------------------------------------
C.EXCLUDE_CEREB = true; % switch whether cerebellum is accounted as subc.
C.atlases     = {'AAL','DBS80','DK84','Tian_S1','Tian_S2','Tian_S3','Tian_S4'};   % fields of `out`
C.tian_scales = {'S1','S2','S3','S4'};
% [OLD] C.dk84_coords = 'aal_tform'; % (version 17/02) | 'dk84_1mm_centre' (version 26/01); used for figure coordinates only

%% ---- Result files (all inside C.P.results_dir) ----------------------------
C.file.results = fullfile(C.P.results_dir, 'subc2rsn_results.mat');        % written by step 3, read by 4-8
C.file.ci      = fullfile(C.P.results_dir, 'subc2rsn_bootstrap_CI.mat');   % written by step 5
C.file.views   = fullfile(C.P.results_dir, 'subc2rsn_views.mat');          % written by step 7a, read by step 7
C.P.out_figures = fullfile(C.P.results_dir, '\Figures');

%% ---- Timeseries smoothing, denoising, filtering (step 1) ----------------------------------------------------
C.DETRND_PreFilt = true;  % switch whether detrend(xRL_dn)
C.ZSCORE_PreFilt = false; % switch whether zscore(detrend(xRL_dn))
C.SAVE_xRL_dn = false;    % switch whether xRL_dn are saved in external disc (or only xRL_dn_f)

C.FWHM=6;      % rfMRI spatial Gaussian smoothing (6mm for T3 images)
C.voxelsize=2; % mm^3
%  Parameters rfMRI butterworth temporal filter
C.TR=0.72;     % Repetition Time (seconds)
C.fnq=1/(2*C.TR);                % Nyquist frequency
C.flp = 0.008;                   % low boundary (highpass) frequency of filter (Hz)
C.fhi = 0.08;                    % high boundary (lowpass)
 % Discarded time points after temporal filter for FC (step 2, cluster)
C.discard_t_filt = 0;

%% ---- Cluster (step 2)
C.OVERLAPPING_VXLS  = true;       % usage for RSN average FC(subc,RSN) of voxels overlapping as schaefer & subcortical for some atlas

%% ---- Bootstrap (step 5) ----------------------------------------------------
C.boot.B = 5000;  C.boot.SEED = 7;  C.boot.ALPHA = 0.05;

%% ---- Figure style ----------------------------------------------------------
C.fig.fontName    = 'Helvetica Neue';
C.fig.export_svg  = true;        % write .svg into C.P.out_figures
C.fig.export_fig  = true;       % also write .fig
C.fig.close_after = false;       % close each figure after exporting it

% pies (steps 6)
C.pie.diam_mm     = 25;          % physical diameter of every pie
C.pie.pct_min     = 5;           % percentage shown only if larger than this
C.pie.fs_pct      = 7;   C.pie.fs_title = 8;   C.pie.fs_row = 9;   C.pie.fs_key = 10;

% voxel figures (step 7)
C.vox.figW_mm = 180;  C.vox.figH_mm = 90;
C.vox.bckg_size = 10; C.vox.roi_size = 15;
C.vox.tickFS = 6;  C.vox.labelFS = 8;  C.vox.titleFS = 8;
C.vox.view_S1_L = [ 90 15];  C.vox.view_S1_R = [-90 15];  C.vox.view_S4 = [90 15];   % used when no saved view exists

% alluvial flows (step 8): options passed to plotSankeyFlowChart_memory_yeo_cm
C.alluvial.options = struct('flow_transparency', 0.3, 'bar_width', 20, 'show_perc', true, 'text_color', 'k', ...
                            'show_cat_labels', true, 'show_layer_labels', true, 'show_legend', false);
C.alluvial.fontSize = 12;
C.alluvial.rsn_order = 1:7; % top-to-bottom order of the networks, any key

C.alluvial.rsn_order_key = struct(); % OR if we want per-key optional order
C.alluvial.rsn_order_key.HIP = [2 7 1 3]; %  e.g. the one tuned by hand:

end
