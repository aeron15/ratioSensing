function plot_figure1B

% Plot the data for figure 1
gal_conc=-9:0.5:2;
glc_conc=-9:0.5:0;
gal_final = [0 2.^[gal_conc]];
glc_final = [0 2.^[glc_conc]];

gal_m = 2.^[gal_conc(1)-1:0.5:gal_conc(end)+0.5];
glc_m = 2.^[glc_conc(1)-1:0.5:glc_conc(end)+0.5];

load('../data/20140701_stitched_areas/output/M_for_replicate_original_supp_material.mat')
load('../data/20140701_stitched_areas/output/M_D_stitched.mat')
%% 20140708 Supplementary different cutoffs

compare_thresholds_fig1B(gal_final,glc_final,gal_m,glc_m,M,D)
compare_thresholds_fig_fitness_schematic(gal_final,glc_final,gal_m,glc_m,M,D)

%% Single responses for the WT strain

% load('output/M_D_stitched.mat')
% 
%predict_responses_renan(gal_final,glc_final,gal_m,glc_m,D,M)

end