function driver_plot_figures()

%% Plot decision fronts for panel E of figure 1
AnalyzeInductionData()
AnalyzeInductionData_A10_H3()

%% Plots mig1? and WT data, also gal80? and WT
%alaso adds the decision fronts to the deletes (but decision front is
%missing)
% figure3
% figure3_YS

%% Generates the WT and the gal2 delete experiment 
%figure4

%% Generates the different models of ratio sensing and how they can be achieved. Fix output data
%figure2

%% Plots the multiplication of the single input functions, I generated a separate version with the stitched data.
% This script also plots the effect of depletion
%predictResponse()

%% Plots the data of the heterozygous deletes
%figureSI4


%% Figure SI3. Plot data at different inoclum sizes and different decision fronts
%figureSI3

%% Compares the different metrics (area, mean, percentage) for an s288c replicate and for the BC187 (A10) strain. For supplementary table.
% Check the version for the stitched data with lower concentrations. Many
% figures of the main paper are here too.
%CompareMetrics


%% %Makes S plot and replicates. Plots the histograms for figure 1 too.
%s282cSurfaceFit

%% Compare Mig1 localization to gal80?
%CompareMig1localiztionToGal80

%% Plot previously assay concentrations of galactose and glucose
PlotDecision_Fronts_Table()
end