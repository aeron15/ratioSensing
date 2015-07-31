function driver_plot_SUPP_figures()

%% Fig. S2. Measured response compared to a multiplicative prediction for S288C
plot_figureS2()

%% Fig. S6. Different metrics for quantifying the response of S288C
plot_figure_S6

plot_figure_S6C_S_plot

%>>>>>>>>>NEED TO REVISE CODE BELOW

plot_FigureS1
%% Plot figure 1C. Plot S for separate replicates. Probably need the one for the data on the figure.
%s288cSurfaceFit
%% Generates the different models of ratio sensing and how they can be achieved. Fix output data
%figure2

%% Plots the multiplication of the single input functions, I generated a separate version with the stitched data.
% This script also plots the effect of depletion
%predictResponse()
%predict_response_multiplication()%from stitched data

%% Plots the data of the heterozygous deletes
%figureSI4

%% Figure SI3. Plot data at different inoclum sizes and different decision fronts
%figureSI3

%% Compares the different metrics (area, mean, percentage) for an s288c replicate and for the BC187 (A10) strain. For supplementary table.
% Check the version for the stitched data with lower concentrations. Many
% figures of the main paper are here too.
%CompareMetrics

%% Compare Mig1 localization to gal80?
%CompareMig1localiztionToGal80

%% Plot previously assay concentrations of galactose and glucose. The lowest concentrations are actually NO  concentrations.
%PlotDecision_Fronts_Table()
end