% pipepline to generate all figures 

addpath('/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/functions/')

% select path where to save figures
Path_figures='/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/figures';

% select path where data are stored
Path_variables='/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/variables';
Behav_path='/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/data/behavior';
MEG_path='/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/data/meg';


%figure 1

Figure_1;

% get significant sensors

get_significant_channels;

% calculate TFR at sensor level;

TFR_sensor_level;

% TFR_source_space

TFR_source_space;

% plot  Figure 2

Figure_2;

%  TFR_source_space for each drug

TFR_source_space_drugs;

% plot Figure 3 and 3.1;
Code_path='/Users/tosot/Downloads/ASSR_newfig/pdf_fig/revisions/final/codes/';


Figure_3;

% ITPC TFR_source_space for PL 40 hz

ITPC_source_space;


% ITPC TFR_source_space for each drug

ITPC_source_space_drugs;

% plot Figure  3-2

Figure_3_2;

