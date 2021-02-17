%% Figure 2
clc; clear all; close all;
% Yearly Met variables
Path_Met_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Historical\Met_Var_Historical_Met_Year.mat';
Path_Met_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp126\Met_Var_ssp126_Met_Year.mat';
Path_Met_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp245\Met_Var_ssp245_Met_Year.mat';
Path_Met_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp370\Met_Var_ssp370_Met_Year.mat';
Path_Met_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp585\Met_Var_ssp585_Met_Year.mat';
Path_Met_Year_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Princeton\Met_Year_Princeton.mat';

Met_Year = cat(2,load(Path_Met_Year_Historical_GCM),...
    load(Path_Met_Year_ScenarioMIP_ssp126_GCM),load(Path_Met_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Met_Year_ScenarioMIP_ssp370_GCM),load(Path_Met_Year_ScenarioMIP_ssp585_GCM),...
    load(Path_Met_Year_Princeton_GCM));

clear Path_Met_Year_Historical_GCM Path_Met_Year_ScenarioMIP_ssp126_GCM Path_Met_Year_ScenarioMIP_ssp245_GCM
clear Path_Met_Year_ScenarioMIP_ssp370_GCM Path_Met_Year_ScenarioMIP_ssp585_GCM Path_Met_Year_Princeton_GCM

Fig2_Left_Plotting(Met_Year)
Fig2_Right_Plotting(Met_Year)