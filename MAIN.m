function MAIN

%MAIN script creates all the figures from the PNAS paper

%Create an output folder if it does not exist
if ~exist('../output_ratio_sensing')
    mkdir('../output_ratio_sensing');
end

if ~exist('../output_ratio_sensing')
    mkdir('../output_ratio_sensing');
end

path_data='../output_ratio_sensing/';


%% Plot main figures of the paper
driver_plot_MAIN_figures()

%% Plot supplementary figures

driver_plot_SUPP_figures()
 
%% generates the uptake heat maps
Transporter_competition