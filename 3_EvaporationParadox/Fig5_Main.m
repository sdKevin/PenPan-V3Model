%% Figure 3
% Princeton 1948-2014, CMIP6 Historical 1850-2014, CMIP6 Scenarios 2015-2100
%% Spatial Distribution of 601B
clc; clear all; close all;
Path_Ensemble_Mean_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\Historical\Epan_Historical_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp126\Epan_ssp126_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp245\Epan_ssp245_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp370\Epan_ssp370_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp585\Epan_ssp585_Ensemble_Mean.mat';
Path_Ensemble_Mean_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\Princeton\Ensemble_Mean.mat';
Path_Fig5_Output = 'Fig5_Output\601B\';

GridEpan_CMIP = cat(2,load(Path_Ensemble_Mean_Princeton_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM),...
    load(Path_Ensemble_Mean_Historical_GCM));

clear Path_Ensemble_Mean_Historical_GCM Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM
clear Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM Path_Ensemble_Mean_Princeton_GCM
% Epan Input (m/s) Output mm/(year2)
Fig5a_Plotting(GridEpan_CMIP , Path_Fig5_Output)

%% Spatial Distribution of D20
clc; clear all; close all;
Path_Ensemble_Mean_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Historical\Epan_Historical_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp126\Epan_ssp126_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp245\Epan_ssp245_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp370\Epan_ssp370_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp585\Epan_ssp585_Ensemble_Mean.mat';
Path_Ensemble_Mean_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Princeton\Ensemble_Mean.mat';
Path_Fig5_Output = 'Fig5_Output\D20\';

GridEpan_CMIP = cat(2,load(Path_Ensemble_Mean_Princeton_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM),...
    load(Path_Ensemble_Mean_Historical_GCM));

clear Path_Ensemble_Mean_Historical_GCM Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM
clear Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM Path_Ensemble_Mean_Princeton_GCM
% Epan Input (m/s) Output mm/(year2)
Fig5a_Plotting(GridEpan_CMIP , Path_Fig5_Output)

%% Spatial Distribution of ClassA
clc; clear all; close all;
Path_Ensemble_Mean_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\Historical\Epan_Historical_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp126\Epan_ssp126_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp245\Epan_ssp245_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp370\Epan_ssp370_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp585\Epan_ssp585_Ensemble_Mean.mat';
Path_Ensemble_Mean_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\Princeton\Ensemble_Mean.mat';
Path_Fig5_Output = 'Fig5_Output\ClassA\';

GridEpan_CMIP = cat(2,load(Path_Ensemble_Mean_Princeton_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM),...
    load(Path_Ensemble_Mean_Historical_GCM));

clear Path_Ensemble_Mean_Historical_GCM Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM
clear Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM Path_Ensemble_Mean_Princeton_GCM
% Epan Input (m/s) Output mm/(year2)
Fig5a_Plotting(GridEpan_CMIP , Path_Fig5_Output)

%% Spatial Distribution of GGI3000
clc; clear all; close all;
Path_Ensemble_Mean_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\Historical\Epan_Historical_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp126\Epan_ssp126_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp245\Epan_ssp245_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp370\Epan_ssp370_Ensemble_Mean.mat';
Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp585\Epan_ssp585_Ensemble_Mean.mat';
Path_Ensemble_Mean_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\Princeton\Ensemble_Mean.mat';
Path_Fig5_Output = 'Fig5_Output\GGI3000\';

GridEpan_CMIP = cat(2,load(Path_Ensemble_Mean_Princeton_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM),...
    load(Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM),load(Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM),...
    load(Path_Ensemble_Mean_Historical_GCM));

clear Path_Ensemble_Mean_Historical_GCM Path_Ensemble_Mean_ScenarioMIP_ssp126_GCM Path_Ensemble_Mean_ScenarioMIP_ssp245_GCM
clear Path_Ensemble_Mean_ScenarioMIP_ssp370_GCM Path_Ensemble_Mean_ScenarioMIP_ssp585_GCM Path_Ensemble_Mean_Princeton_GCM
% Epan Input (m/s) Output mm/(year2)
Fig5a_Plotting(GridEpan_CMIP , Path_Fig5_Output)