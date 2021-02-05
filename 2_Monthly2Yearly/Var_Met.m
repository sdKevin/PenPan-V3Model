clc; clear all; close all;

%% Setting the input/output paths
% CMIP6 Historical Meteorological Data
InputPath{1} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\Historical\Met_Var_Historical_';
OutputPath{1} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Historical\Met_Var_Historical_';
% CMIP6 ScenarioMIP ssp126 Meteorological Data
InputPath{2} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\ScenarioMIP_ssp126\Met_Var_ssp126_';
OutputPath{2} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp126\Met_Var_ssp126_';
% CMIP6 ScenarioMIP ssp245 Meteorological Data
InputPath{3} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\ScenarioMIP_ssp245\Met_Var_ssp245_';
OutputPath{3} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp245\Met_Var_ssp245_';
% CMIP6 ScenarioMIP ssp370 Meteorological Data
InputPath{4} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\ScenarioMIP_ssp370\Met_Var_ssp370_';
OutputPath{4} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp370\Met_Var_ssp370_';
% CMIP6 ScenarioMIP ssp585 Meteorological Data
InputPath{5} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\ScenarioMIP_ssp585\Met_Var_ssp585_';
OutputPath{5} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\ScenarioMIP_ssp585\Met_Var_ssp585_';

%% (1) Integrate Monthly CMIP6 data to yearly data
for i_Path = 1 : length(InputPath)
    if i_Path == 4
        % Name of Global Climate Model
        GCM_Ensemble = {'ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-CanOE',...
            'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg',...
            'FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','INM-CM4-8',...
            'INM-CM5-0','IPSL-CM6A-LR','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
            'MRI-ESM2-0','NorESM2-MM','UKESM1-0-LL'};
    else
        % Name of Global Climate Model, since HadGEM3-GC31-LL model does
        % not have ssp370
        GCM_Ensemble = {'ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-CanOE',...
            'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg',...
            'FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8',...
            'INM-CM5-0','IPSL-CM6A-LR','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
            'MRI-ESM2-0','NorESM2-MM','UKESM1-0-LL'};
    end
    for i_GCM = 1 : length(GCM_Ensemble)
        %% (1.1) 1850-2014 Historical Year Series
        GCM = GCM_Ensemble{i_GCM};
        load(strcat(InputPath{i_Path} , GCM , '.mat'));
        % Monthly Grid to Yearly Grid : GridYear
        iii = 1;
        for ii = 1 : 12 : size(Met_Var.Sg , 3)
            % Sg
            A = Met_Var.Sg(:,:,ii:ii+11);
            GridYear.Sg(:,:,iii) = nanmean(A,3);
            % Ta
            B = Met_Var.Ra(:,:,ii:ii+11);
            GridYear.Ra(:,:,iii) = nanmean(B,3);
            % VPD
            C = Met_Var.Li(:,:,ii:ii+11);
            GridYear.Li(:,:,iii) = nanmean(C,3);
            % CO2
            D = Met_Var.U10(:,:,ii:ii+11);
            GridYear.U10(:,:,iii) = nanmean(D,3);
            % U2
            E = Met_Var.Ta(:,:,ii:ii+11);
            GridYear.Ta(:,:,iii) = nanmean(E,3);
            % pr
            F = Met_Var.Sh(:,:,ii:ii+11);
            GridYear.Sh(:,:,iii) = nanmean(F,3);
            clear A B C D E F
            iii = iii+1;
        end
        clear ii iii;
        % Yearly Grid to Yearly series : Met_Year
        for ii = 1:size(GridYear.Sg,3)
            A = GridYear.Sg(:,:,ii); Met_Year.Sg(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.Ra(:,:,ii); Met_Year.Ra(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.Li(:,:,ii); Met_Year.Li(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.U10(:,:,ii); Met_Year.U10(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.Ta(:,:,ii); Met_Year.Ta(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.Sh(:,:,ii); Met_Year.Sh(i_GCM , ii) = nanmean(A(:)); clear A;
        end
        clear ii Met_Var
        Met_Var = GridYear; clear GridYear;
        % Save Epan to calculate Ensemble Mean
        All_Met_Var.Sg(:,:,:,i_GCM) = Met_Var.Sg;
        All_Met_Var.Ra(:,:,:,i_GCM) = Met_Var.Ra;
        All_Met_Var.Li(:,:,:,i_GCM) = Met_Var.Li;
        All_Met_Var.U10(:,:,:,i_GCM) = Met_Var.U10;
        All_Met_Var.Ta(:,:,:,i_GCM) = Met_Var.Ta;
        All_Met_Var.Sh(:,:,:,i_GCM) = Met_Var.Sh;
        %% (1.2) Output GridYear from 1850-2100
        save(strcat(OutputPath{i_Path} , GCM) , 'Met_Var');
        clear Met_Var GCM
    end
    clear GCM i_GCM
    save(strcat(OutputPath{i_Path} , 'Met_Year') , 'Met_Year');
    
    Ensemble_Mean_Met.Sg = nanmean(All_Met_Var.Sg,4);
    Ensemble_Mean_Met.Ra = nanmean(All_Met_Var.Ra,4);
    Ensemble_Mean_Met.Li = nanmean(All_Met_Var.Li,4);
    Ensemble_Mean_Met.U10 = nanmean(All_Met_Var.U10,4);
    Ensemble_Mean_Met.Ta = nanmean(All_Met_Var.Ta,4);
    Ensemble_Mean_Met.Sh = nanmean(All_Met_Var.Sh,4);
    
    save(strcat(OutputPath{i_Path} , 'Ensemble_Mean') , 'Ensemble_Mean_Met');
    clear Met_Year All_Met_Var
end

%% (2) Integrate Monthly Princeton data to yearly data
clc; clear all; close all;
% CMIP6 Princeton Meteorological Data
InputPath_Princeton = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\Princeton\Met_Var_Princeton.mat';
OutputPath_Princeton = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_Met\Princeton\';

load(InputPath_Princeton);
% Monthly Grid to Yearly Grid : GridYear
iii = 1;
for ii = 1 : 12 : size(Met_Var.Sg , 3)
    % Sg
    A = Met_Var.Sg(:,:,ii:ii+11);
    GridYear.Sg(:,:,iii) = nanmean(A,3);
    % Ta
    B = Met_Var.Ra(:,:,ii:ii+11);
    GridYear.Ra(:,:,iii) = nanmean(B,3);
    % VPD
    C = Met_Var.Li(:,:,ii:ii+11);
    GridYear.Li(:,:,iii) = nanmean(C,3);
    % CO2
    D = Met_Var.U10(:,:,ii:ii+11);
    GridYear.U10(:,:,iii) = nanmean(D,3);
    % U2
    E = Met_Var.Ta(:,:,ii:ii+11);
    GridYear.Ta(:,:,iii) = nanmean(E,3);
    % pr
    F = Met_Var.Sh(:,:,ii:ii+11);
    GridYear.Sh(:,:,iii) = nanmean(F,3);
    clear A B C D E F
    iii = iii+1;
end
clear ii iii;
% Yearly Grid to Yearly series : Met_Year
for ii = 1:size(GridYear.Sg,3)
    A = GridYear.Sg(:,:,ii); Met_Year.Sg(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.Ra(:,:,ii); Met_Year.Ra(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.Li(:,:,ii); Met_Year.Li(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.U10(:,:,ii); Met_Year.U10(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.Ta(:,:,ii); Met_Year.Ta(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.Sh(:,:,ii); Met_Year.Sh(i_GCM , ii) = nanmean(A(:)); clear A;
end
clear ii Met_Var
Met_Var = GridYear; clear GridYear;
Ensemble_Mean_Met = Met_Var;
%% (1.2) Output GridYear from 1850-2100
save(strcat(OutputPath_Princeton , 'Met_Var_Princeton') , 'Met_Var');
clear Met_Var
save(strcat(OutputPath_Princeton , 'Met_Year_Princeton') , 'Met_Year');
clear Met_Year
save(strcat(OutputPath_Princeton , 'Ensemble_Mean') , 'Ensemble_Mean_Met');
clear Ensemble_Mean_Met