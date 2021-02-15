%% Figure 1
clc; clear all; close all;
% Yearly Met variables
Path_Epan_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\Historical\Epan_Historical_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp126\Epan_ssp126_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp245\Epan_ssp245_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp370\Epan_ssp370_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp585\Epan_ssp585_Epan_Year.mat';

Epan_Year = cat(2,load(Path_Epan_Year_Historical_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp126_GCM),load(Path_Epan_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp370_GCM),load(Path_Epan_Year_ScenarioMIP_ssp585_GCM));

clear Path_Epan_Year_Historical_GCM Path_Epan_Year_ScenarioMIP_ssp126_GCM Path_Epan_Year_ScenarioMIP_ssp245_GCM
clear Path_Epan_Year_ScenarioMIP_ssp370_GCM Path_Epan_Year_ScenarioMIP_ssp585_GCM Path_Epan_Year_Princeton_GCM

% Input unit: m/s
Fig1_Left_Plotting(Epan_Year)
Fig1_Right_Plotting(Epan_Year)