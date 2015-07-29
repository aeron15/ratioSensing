function driver_plot_MAIN_figures()

%% Plot figure 1B. Double gradient and decision front.s282cSurfaceFit also plots histograms figure 1.
plot_fig1B()

%% Plot figure 1C.Fraction of inducing cells as a function of the ratio of galactose and glucose concentrations.
plot_figure1C()

%% Plot schematic of figure 4A
plot_fig4A()

%% Plot figure 1C. Plot S for separate replicates. Probably need the one for the data on the figure.
%s288cSurfaceFit

%% Plot decision fronts of BC187, YJM978 and S288C for panel E of figure 1
%AnalyzeInductionData()
%AnalyzeInductionData_A10_H3()

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