function driver_plot_MAIN_figures()
%PLOTS all figures of the main paper PNAS

%% Plot figure 1B. Double gradient and decision front.s282cSurfaceFit also plots histograms figure 1.
plot_fig1B
plot_fig1B_histograms
plot_histograms_cellasic_fig1B

%% Plot figure 1C.Fraction of inducing cells as a function of the ratio of galactose and glucose concentrations.
plot_figure1C

%% Plot figure 1D()
%Comparison of models of signal integration (SI Appendix) by threshold sensing (Upper) and ratio sensing (Lower), displayed as in B and C.
plot_1D

%% Plot decision fronts of BC187, YJM978 and S288C for panel E of figure 1.
%Decision fronts, calculated as in B, for three strains of S. cerevisiae.
plot_figure1E

%% Plots gal80 delete
plot_Figure2C

%% Plot mig1 delete
plot_Figure2D

%% Plot gal80 delete mig1 delete
plot_Figure2E

%% Generates the WT and the gal2 delete experiment for figure 3
%Deletion of GAL2 does not eliminate ratio sensing. Black and red lines are the decision front
plot_Figure3

%% Plot schematic of figure 4A
plot_fig4A

%% Plot modes of signal integration
plot_Figure4B

end