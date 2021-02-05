clc; clear all; close all;

%% Setting the input/output paths
% CMIP6 Historical Meteorological Data
InputPath{1} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\Historical\Epan_Historical_';
OutputPath{1} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Historical\Epan_Historical_';
% CMIP6 ScenarioMIP ssp126 Meteorological Data
InputPath{2} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\ScenarioMIP_ssp126\Epan_ssp126_';
OutputPath{2} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp126\Epan_ssp126_';
% CMIP6 ScenarioMIP ssp245 Meteorological Data
InputPath{3} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\ScenarioMIP_ssp245\Epan_ssp245_';
OutputPath{3} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp245\Epan_ssp245_';
% CMIP6 ScenarioMIP ssp370 Meteorological Data
InputPath{4} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\ScenarioMIP_ssp370\Epan_ssp370_';
OutputPath{4} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp370\Epan_ssp370_';
% CMIP6 ScenarioMIP ssp585 Meteorological Data
InputPath{5} = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\ScenarioMIP_ssp585\Epan_ssp585_';
OutputPath{5} = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\ScenarioMIP_ssp585\Epan_ssp585_';

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
        % Monthly Grid to Yearly Grid : GridYear [m/s]
        iii = 1;
        for ii = 1 : 12 : size(Epan.E_pan,3)
            A = Epan.E_pan(:,:,ii:ii+11);
            GridYear.E_pan(:,:,iii) = nanmean(A,3);
            B = Epan.E_pan_R(:,:,ii:ii+11);
            GridYear.E_pan_R(:,:,iii) = nanmean(B,3);
            C = Epan.E_pan_A(:,:,ii:ii+11);
            GridYear.E_pan_A(:,:,iii) = nanmean(C,3);
            D = Epan.E_pan_R_w(:,:,ii:ii+11);
            GridYear.E_pan_R_w(:,:,iii) = nanmean(D,3);
            E = Epan.E_pan_R_wall(:,:,ii:ii+11);
            GridYear.E_pan_R_wall(:,:,iii) = nanmean(E,3);
            F = Epan.E_pan_R_rim(:,:,ii:ii+11);
            GridYear.E_pan_R_rim(:,:,iii) = nanmean(F,3);
            AA = Epan.E_pan_R_rim_pie(:,:,ii:ii+11);
            GridYear.E_pan_R_rim_pie(:,:,iii) = nanmean(AA,3);
            BB = Epan.E_pan_R_bot(:,:,ii:ii+11);
            GridYear.E_pan_R_bot(:,:,iii) = nanmean(BB,3);
            clear A B C D E F AA BB CC DD EE
            iii = iii+1;
        end
        clear ii iii;
        % Yearly Grid to Yearly series : Epan_Year
        for ii = 1 : size(GridYear.E_pan , 3)
            A = GridYear.E_pan(:,:,ii); Epan_Year.E_pan(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R(:,:,ii); Epan_Year.E_pan_R(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_A(:,:,ii); Epan_Year.E_pan_A(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R_w(:,:,ii); Epan_Year.E_pan_R_w(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R_wall(:,:,ii); Epan_Year.E_pan_R_wall(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R_rim(:,:,ii); Epan_Year.E_pan_R_rim(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R_rim_pie(:,:,ii); Epan_Year.E_pan_R_rim_pie(i_GCM , ii) = nanmean(A(:)); clear A;
            A = GridYear.E_pan_R_bot(:,:,ii); Epan_Year.E_pan_R_bot(i_GCM , ii) = nanmean(A(:)); clear A;
        end
        clear ii Epan
        Epan = GridYear; clear GridYear;
        % Save Epan to calculate Ensemble Mean
        All_Epan.E_pan(:,:,:,i_GCM) = Epan.E_pan;
        All_Epan.E_pan_R(:,:,:,i_GCM) = Epan.E_pan_R;
        All_Epan.E_pan_A(:,:,:,i_GCM) = Epan.E_pan_A;
        All_Epan.E_pan_R_w(:,:,:,i_GCM) = Epan.E_pan_R_w;
        All_Epan.E_pan_R_wall(:,:,:,i_GCM) = Epan.E_pan_R_wall;
        All_Epan.E_pan_R_rim(:,:,:,i_GCM) = Epan.E_pan_R_rim;
        All_Epan.E_pan_R_rim_pie(:,:,:,i_GCM) = Epan.E_pan_R_rim_pie;
        All_Epan.E_pan_R_bot(:,:,:,i_GCM) = Epan.E_pan_R_bot;
        %% (1.2) Output GridYear from 1850-2100
        save(strcat(OutputPath{i_Path} , GCM) , 'Epan');
        clear Epan GCM
    end
    clear i_GCM
    save(strcat(OutputPath{i_Path} , 'Epan_Year') , 'Epan_Year');
    
    Ensemble_Mean_Epan.E_pan = nanmean(All_Epan.E_pan,4);
    Ensemble_Mean_Epan.E_pan_R = nanmean(All_Epan.E_pan_R,4);
    Ensemble_Mean_Epan.E_pan_A = nanmean(All_Epan.E_pan_A,4);
    Ensemble_Mean_Epan.E_pan_R_w = nanmean(All_Epan.E_pan_R_w,4);
    Ensemble_Mean_Epan.E_pan_R_wall = nanmean(All_Epan.E_pan_R_wall,4);
    Ensemble_Mean_Epan.E_pan_R_rim = nanmean(All_Epan.E_pan_R_rim,4);
    Ensemble_Mean_Epan.E_pan_R_rim_pie = nanmean(All_Epan.E_pan_R_rim_pie,4);
    Ensemble_Mean_Epan.E_pan_R_bot = nanmean(All_Epan.E_pan_R_bot,4);
    
    save(strcat(OutputPath{i_Path} , 'Ensemble_Mean') , 'Ensemble_Mean_Epan');
    clear Epan_Year All_Epan
end

%% (2) Integrate Monthly Princeton data to yearly data
clc; clear all; close all;
% Princeton Meteorological Data
InputPath_Princeton = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan\Princeton\Epan_Princeton.mat';
OutputPath_Princeton = 'E:\PenPanV3\VariableStorage\YearlyVar\Var_D20Epan\Princeton\';

load(InputPath_Princeton);
% Monthly Grid to Yearly Grid : GridYear
iii = 1;
for ii = 1 : 12 : size(Epan.E_pan,3)
    A = Epan.E_pan(:,:,ii:ii+11);
    GridYear.E_pan(:,:,iii) = nanmean(A,3);
    B = Epan.E_pan_R(:,:,ii:ii+11);
    GridYear.E_pan_R(:,:,iii) = nanmean(B,3);
    C = Epan.E_pan_A(:,:,ii:ii+11);
    GridYear.E_pan_A(:,:,iii) = nanmean(C,3);
    D = Epan.E_pan_R_w(:,:,ii:ii+11);
    GridYear.E_pan_R_w(:,:,iii) = nanmean(D,3);
    E = Epan.E_pan_R_wall(:,:,ii:ii+11);
    GridYear.E_pan_R_wall(:,:,iii) = nanmean(E,3);
    F = Epan.E_pan_R_rim(:,:,ii:ii+11);
    GridYear.E_pan_R_rim(:,:,iii) = nanmean(F,3);
    AA = Epan.E_pan_R_rim_pie(:,:,ii:ii+11);
    GridYear.E_pan_R_rim_pie(:,:,iii) = nanmean(AA,3);
    BB = Epan.E_pan_R_bot(:,:,ii:ii+11);
    GridYear.E_pan_R_bot(:,:,iii) = nanmean(BB,3);
    clear A B C D E F AA BB CC DD EE
    iii = iii+1;
end
clear ii iii;
% Yearly Grid to Yearly series : Epan_Year
for ii = 1 : size(GridYear.E_pan , 3)
    A = GridYear.E_pan(:,:,ii); Epan_Year.E_pan(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R(:,:,ii); Epan_Year.E_pan_R(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_A(:,:,ii); Epan_Year.E_pan_A(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R_w(:,:,ii); Epan_Year.E_pan_R_w(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R_wall(:,:,ii); Epan_Year.E_pan_R_wall(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R_rim(:,:,ii); Epan_Year.E_pan_R_rim(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R_rim_pie(:,:,ii); Epan_Year.E_pan_R_rim_pie(i_GCM , ii) = nanmean(A(:)); clear A;
    A = GridYear.E_pan_R_bot(:,:,ii); Epan_Year.E_pan_R_bot(i_GCM , ii) = nanmean(A(:)); clear A;
end
clear ii Epan
Epan = GridYear; clear GridYear;
Ensemble_Mean_Epan = Epan;
%% (1.2) Output GridYear from 1850-2100
save(strcat(OutputPath_Princeton , 'Epan_Princeton') , 'Epan');
clear Epan
save(strcat(OutputPath_Princeton , 'Epan_Year_Princeton') , 'Epan_Year');
clear Epan_Year
save(strcat(OutputPath_Princeton , 'Ensemble_Mean') , 'Ensemble_Mean_Epan');
clear Ensemble_Mean_Epan