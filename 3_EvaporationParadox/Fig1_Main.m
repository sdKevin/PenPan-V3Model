%% Figure 1 601B pan
clc; clear all; close all;
% Yearly Epan of 601B pan
Path_Epan_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\Historical\Epan_Historical_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp126\Epan_ssp126_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp245\Epan_ssp245_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp370\Epan_ssp370_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\ScenarioMIP_ssp585\Epan_ssp585_Epan_Year.mat';
Path_Epan_Year_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_601BEpan\Princeton\Epan_Year_Princeton.mat';

Epan_Year = cat(2,load(Path_Epan_Year_Historical_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp126_GCM),load(Path_Epan_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp370_GCM),load(Path_Epan_Year_ScenarioMIP_ssp585_GCM),...
    load(Path_Epan_Year_Princeton_GCM));

clear Path_Epan_Year_Historical_GCM Path_Epan_Year_ScenarioMIP_ssp126_GCM Path_Epan_Year_ScenarioMIP_ssp245_GCM
clear Path_Epan_Year_ScenarioMIP_ssp370_GCM Path_Epan_Year_ScenarioMIP_ssp585_GCM Path_Epan_Year_Princeton_GCM

% Input unit: m/s
YlimRange_Epan = [-50,300]; YTickRange_Epan = [-50:50:300]; YTickLabel_Epan = {'';'0';'';'100';'';'200';'';'300'};
YlimRange_Epan_R = [-40,160]; YTickRange_Epan_R = [-40:40:160]; YTickLabel_Epan_R = {'-40';'0';'40';'80';'120';'160'};
YlimRange_Epan_A = [-20,100]; YTickRange_Epan_A = [-20:20:100]; YTickLabel_Epan_A = {'-20';'0';'20';'40';'60';'80';'100'};

Fig1_Left_Plotting(Epan_Year,YlimRange_Epan,YTickRange_Epan,YTickLabel_Epan,...
    YlimRange_Epan_R,YTickRange_Epan_R,YTickLabel_Epan_R,YlimRange_Epan_A,YTickRange_Epan_A,YTickLabel_Epan_A)
Fig1_Right_Plotting(Epan_Year,YlimRange_Epan,YlimRange_Epan_R,YlimRange_Epan_A)

%% Figure 1 D20 pan
clc; clear all; close all;
% Yearly Epan of D20 pan
Path_Epan_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Historical\Epan_Historical_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp126\Epan_ssp126_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp245\Epan_ssp245_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp370\Epan_ssp370_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp585\Epan_ssp585_Epan_Year.mat';
Path_Epan_Year_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Princeton\Epan_Year_Princeton.mat';

Epan_Year = cat(2,load(Path_Epan_Year_Historical_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp126_GCM),load(Path_Epan_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp370_GCM),load(Path_Epan_Year_ScenarioMIP_ssp585_GCM),...
    load(Path_Epan_Year_Princeton_GCM));

clear Path_Epan_Year_Historical_GCM Path_Epan_Year_ScenarioMIP_ssp126_GCM Path_Epan_Year_ScenarioMIP_ssp245_GCM
clear Path_Epan_Year_ScenarioMIP_ssp370_GCM Path_Epan_Year_ScenarioMIP_ssp585_GCM Path_Epan_Year_Princeton_GCM

% Input unit: m/s
YlimRange_Epan = [-100,800]; YTickRange_Epan = [-100:100:800]; YTickLabel_Epan = {'';'0';'';'200';'';'400';'';'600';'';'800'};
YlimRange_Epan_R = [-50,300]; YTickRange_Epan_R = [-50:50:300]; YTickLabel_Epan_R = {'';'0';'';'100';'';'200';'';'300'};
YlimRange_Epan_A = [-50,500]; YTickRange_Epan_A = [-50:50:500]; YTickLabel_Epan_A = {'';'0';'';'100';'';'200';'';'300';'';'400';'';'500'};

Fig1_Left_Plotting(Epan_Year,YlimRange_Epan,YTickRange_Epan,YTickLabel_Epan,...
    YlimRange_Epan_R,YTickRange_Epan_R,YTickLabel_Epan_R,YlimRange_Epan_A,YTickRange_Epan_A,YTickLabel_Epan_A)
Fig1_Right_Plotting(Epan_Year,YlimRange_Epan,YlimRange_Epan_R,YlimRange_Epan_A)

%% Figure 1 ClassA pan
clc; clear all; close all;
% Yearly Epan of ClassA pan
Path_Epan_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\Historical\Epan_Historical_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp126\Epan_ssp126_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp245\Epan_ssp245_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp370\Epan_ssp370_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\ScenarioMIP_ssp585\Epan_ssp585_Epan_Year.mat';
Path_Epan_Year_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_ClassAEpan\Princeton\Epan_Year_Princeton.mat';

Epan_Year = cat(2,load(Path_Epan_Year_Historical_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp126_GCM),load(Path_Epan_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp370_GCM),load(Path_Epan_Year_ScenarioMIP_ssp585_GCM),...
    load(Path_Epan_Year_Princeton_GCM));

clear Path_Epan_Year_Historical_GCM Path_Epan_Year_ScenarioMIP_ssp126_GCM Path_Epan_Year_ScenarioMIP_ssp245_GCM
clear Path_Epan_Year_ScenarioMIP_ssp370_GCM Path_Epan_Year_ScenarioMIP_ssp585_GCM Path_Epan_Year_Princeton_GCM

% Input unit: m/s
YlimRange_Epan = [-100,500]; YTickRange_Epan = [-100:100:500]; YTickLabel_Epan = {'-100';'0';'100';'200';'300';'400';'500';};
YlimRange_Epan_R = [-50,200]; YTickRange_Epan_R = [-50:50:200]; YTickLabel_Epan_R = {'-50';'0';'50';'100';'150';'200';};
YlimRange_Epan_A = [-50,300]; YTickRange_Epan_A = [-50:50:300]; YTickLabel_Epan_A = {'';'0';'';'100';'';'200';'';'300';};

Fig1_Left_Plotting(Epan_Year,YlimRange_Epan,YTickRange_Epan,YTickLabel_Epan,...
    YlimRange_Epan_R,YTickRange_Epan_R,YTickLabel_Epan_R,YlimRange_Epan_A,YTickRange_Epan_A,YTickLabel_Epan_A)
Fig1_Right_Plotting(Epan_Year,YlimRange_Epan,YlimRange_Epan_R,YlimRange_Epan_A)

%% Figure 1 GGI3000 pan
clc; clear all; close all;
% Yearly Epan of GGI3000 pan
Path_Epan_Year_Historical_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\Historical\Epan_Historical_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp126_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp126\Epan_ssp126_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp245_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp245\Epan_ssp245_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp370_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp370\Epan_ssp370_Epan_Year.mat';
Path_Epan_Year_ScenarioMIP_ssp585_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\ScenarioMIP_ssp585\Epan_ssp585_Epan_Year.mat';
Path_Epan_Year_Princeton_GCM = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_GGI3000Epan\Princeton\Epan_Year_Princeton.mat';

Epan_Year = cat(2,load(Path_Epan_Year_Historical_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp126_GCM),load(Path_Epan_Year_ScenarioMIP_ssp245_GCM),...
    load(Path_Epan_Year_ScenarioMIP_ssp370_GCM),load(Path_Epan_Year_ScenarioMIP_ssp585_GCM),...
    load(Path_Epan_Year_Princeton_GCM));

clear Path_Epan_Year_Historical_GCM Path_Epan_Year_ScenarioMIP_ssp126_GCM Path_Epan_Year_ScenarioMIP_ssp245_GCM
clear Path_Epan_Year_ScenarioMIP_ssp370_GCM Path_Epan_Year_ScenarioMIP_ssp585_GCM Path_Epan_Year_Princeton_GCM

% Input unit: m/s
YlimRange_Epan = [-50,300]; YTickRange_Epan = [-50:50:300]; YTickLabel_Epan = {'';'0';'';'100';'';'200';'';'300'};
YlimRange_Epan_R = [-30,200]; YTickRange_Epan_R = [-30:30:200]; YTickLabel_Epan_R = {'';'0';'';'60';'';'120';'';'180'};
YlimRange_Epan_A = [-20,100]; YTickRange_Epan_A = [-20:20:100]; YTickLabel_Epan_A = {'';'0';'';'40';'';'80';''};

Fig1_Left_Plotting(Epan_Year,YlimRange_Epan,YTickRange_Epan,YTickLabel_Epan,...
    YlimRange_Epan_R,YTickRange_Epan_R,YTickLabel_Epan_R,YlimRange_Epan_A,YTickRange_Epan_A,YTickLabel_Epan_A)
Fig1_Right_Plotting(Epan_Year,YlimRange_Epan,YlimRange_Epan_R,YlimRange_Epan_A)