%% Figure 4
%% Spatial change of Met Variables linear trend
% Princeton 1948-2014, CMIP6 Historical 1850-2014, CMIP6 Scenarios 2015-2100
clc; clear all; close all;
Path_Ensemble_Mean_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Historical\Met_Var_Historical_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp126\Met_Var_ssp126_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp245\Met_Var_ssp245_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp370\Met_Var_ssp370_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp585\Met_Var_ssp585_Ensemble_Mean.mat';
Path_Ensemble_Mean_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Princeton\Ensemble_Mean.mat';
Path_Fig4_Output = 'Fig4_Output\';

GridMet_CMIP = cat(2,load(Path_Ensemble_Mean_Princeton_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM),...
    load(Path_Ensemble_Mean_Historical_GCM));

clear Path_Ensemble_Mean_Historical_GCM Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM
clear Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM Path_Ensemble_Mean_Princeton_GCM
% Input Sg(W/m2),Li(W/m2),U10(m/s),Ta(K),Sh(-) Output Unit/(year)
Fig4a_Plotting(GridMet_CMIP , Path_Fig4_Output)

