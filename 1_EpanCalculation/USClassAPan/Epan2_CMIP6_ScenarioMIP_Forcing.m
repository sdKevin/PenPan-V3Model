clc; clear all; close all;

%% (1) Setting the input/output paths
% CMIP6 Ensemble Meteorological Data
InputPath_CMIP6_Ensemble = 'D:\CMIP6\ProcessData\Ensemble_Met';
% Princeton-GMFD Data
InputPath_Princeton = 'D:\CMIP6\ProcessData\Princeton\monthly';
% Save Pan evaporation (Epan) Data
OutputPath_Epan = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_ClassAEpan';

%% (2) ScenarioMIP Experiment
GCM_Ensemble = {'ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-CanOE',...
    'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg',...
    'FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8',...
    'INM-CM5-0','IPSL-CM6A-LR','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
    'MRI-ESM2-0','NorESM2-MM','UKESM1-0-LL'}; % Name of Global Climate Model
ssp = input('Choose a Scenario (ssp126/ssp245/ssp370/ssp585) : ','s');
for i_GCM = 1 : length(GCM_Ensemble)
    if mean(ssp=='ssp370')==1 && i_GCM==16
        continue;
    end
    %% (2.1) Generate bias correction coefficient
    %% (2.1.1) load CMIP6 Historical Data
    GCM = GCM_Ensemble{i_GCM}
    load(strcat(InputPath_CMIP6_Ensemble , '\Historical\' , GCM , '.mat'));
    
    %% (2.1.2) Interpolating Forcing Data to Uniform Resolution i.e., 0.5deg
    % Load Global 0.5 Degree Coordinate Data from Princeton-GMFD Data
    load LandInfo_05deg;
    % Bilinear Interpolation
    R1.lat = lat_05deg; R1.lon = lon_05deg;
    for ii = 1 : size(r1.huss,3)
        R1.huss(:,:,ii) = interp2(r1.lat , r1.lon , r1.huss(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.ps(:,:,ii) = interp2(r1.lat , r1.lon , r1.ps(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.rlds(:,:,ii) = interp2(r1.lat , r1.lon , r1.rlds(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.rsds(:,:,ii) = interp2(r1.lat , r1.lon , r1.rsds(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.rsdt(:,:,ii) = interp2(r1.lat , r1.lon , r1.rsdt(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.sfcWind(:,:,ii) = interp2(r1.lat , r1.lon , r1.sfcWind(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.tas(:,:,ii) = interp2(r1.lat , r1.lon , r1.tas(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
    end
    clear r1 ii
    r1 = R1; clear R1
    
    %% (2.1.3) Generate Bias Correction Coefficient
    load([InputPath_Princeton , '\huss.mat']); huss(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.huss = nanmean(huss,3) ./ nanmean(r1.huss(:,:,1177:end),3); C.huss(C.huss>10)=10;clear huss;
    load([InputPath_Princeton , '\ps.mat']); ps(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.ps = nanmean(ps,3) ./ nanmean(r1.ps(:,:,1177:end),3); C.ps(C.ps>10)=10;clear ps;
    load([InputPath_Princeton , '\rlds.mat']); rlds(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.rlds = nanmean(rlds,3) ./ nanmean(r1.rlds(:,:,1177:end),3); C.rlds(C.rlds>10)=10;clear rlds;
    load([InputPath_Princeton , '\rsds.mat']); rsds(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.rsds = nanmean(rsds,3) ./ nanmean(r1.rsds(:,:,1177:end),3); C.rsds(C.rsds>10)=10;clear rsds;
    load([InputPath_Princeton , '\sfcWind.mat']); sfcWind(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.sfcWind = nanmean(sfcWind,3) ./ nanmean(r1.sfcWind(:,:,1177:end),3); C.sfcWind(C.sfcWind>10)=10;clear sfcWind;
    load([InputPath_Princeton , '\tas.mat']); tas(:,:,805:828)=[];%1948-2016 to 1948-2014
    C.tas = nanmean(tas,3) - nanmean(r1.tas(:,:,1177:end),3);clear tas;
    clear r1
    
    %% (2.2) Load ssp Data
    load(strcat(InputPath_CMIP6_Ensemble , '\ScenarioMIP\' , ssp , '\',GCM,'.mat'))
    
    %% (2.3) Interpolating Forcing Data to Uniform Resolution i.e., 0.5deg
    % Load Global 0.5 Degree Coordinate Data from Princeton-GMFD Data
    load LandInfo_05deg
    R1.lat = lat_05deg;
    R1.lon = lon_05deg;
    for ii = 1 : size(r1.huss,3)
        R1.huss(:,:,ii) = interp2(r1.lat,r1.lon,r1.huss(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.ps(:,:,ii) = interp2(r1.lat,r1.lon,r1.ps(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.rlds(:,:,ii) = interp2(r1.lat,r1.lon,r1.rlds(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.rsds(:,:,ii) = interp2(r1.lat,r1.lon,r1.rsds(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.rsdt(:,:,ii) = interp2(r1.lat,r1.lon,r1.rsdt(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.sfcWind(:,:,ii) = interp2(r1.lat,r1.lon,r1.sfcWind(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
        R1.tas(:,:,ii) = interp2(r1.lat,r1.lon,r1.tas(:,:,ii),lat_05deg,lon_05deg) .* landmask_05deg;
    end
    clear r1 ii
    r1 = R1; clear R1
    
    %% (2.4) Bias Correction : Correct CMIP6 ScenarioMIP Data
    for ii = 1:size(r1.huss,3)
        r1.huss(:,:,ii) = r1.huss(:,:,ii) .* C.huss;
        r1.ps(:,:,ii) = r1.ps(:,:,ii) .* C.ps;
        r1.rlds(:,:,ii) = r1.rlds(:,:,ii) .* C.rlds;
        r1.rsds(:,:,ii) = r1.rsds(:,:,ii) .* C.rsds;
        r1.rsdt(:,:,ii) = r1.rsdt(:,:,ii) .* C.rsds; % PenPan-V3 uses rsds/rsdt in calculation, thus we use C.rsds to correct rsdt to guaratee the ratio did not change
        r1.sfcWind(:,:,ii) = r1.sfcWind(:,:,ii) .* C.sfcWind;
        r1.tas(:,:,ii) = r1.tas(:,:,ii) + C.tas;
    end
    clear ii
    
    %% (2.5) Calculation pan evaporation (Epan) (m/s)
    %% Pan Parameters
    % L is the diameter of Class A
    pan_pars.D = 1.21; % [m]
    pan_pars.L = pan_pars.D; % [m]
    % he is the height of rim
    pan_pars.he = 0.055; % [m]
    % hw is the height of water level
    pan_pars.hw = 0.2; % [m]
    % Beta is the ratio of heat to mass transfer coefficients of the pan
    pan_pars.Beta = 2 +  pi*pan_pars.D*pan_pars.hw./(0.25*pi*pan_pars.D^2) +  2*pi*pan_pars.D*pan_pars.he./(0.25*pi*pan_pars.D^2);
    % C is the correction factor to account for the shading effect of the bird guard
    pan_pars.C = 1.07; %Class A C=1.07; D20 C=1; 601B C=1;
    % e_gnd is the emissivity of ground
    pan_pars.e_gnd = 0.90;
    % e_wall is the emissivity of water
    pan_pars.e_w = 0.89;
    % e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
    pan_pars.e_wall = 0.82; %Both D20 and US Class A pans are stainless steel coated with ZnO
    % N is the refractive index
    pan_pars.N = 1.33; % Class A: 1.33
    % K is the extinction coeeficient
    pan_pars.K = 0; % Class A: 0
    % The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
    pan_pars.alpha_0_wall = 0.36;
    % alpha_gnd
    pan_pars.alpha_gnd = 0.2;
    
    Epan = PenPan_V3_ClassA(pan_pars , lat_05deg , elevation_05deg ,...
        r1.rsds , r1.rsdt , r1.rlds , r1.sfcWind , r1.tas , r1.huss , r1.ps); % Epan(m/s)
    
    %% (2.6) Save the result
    save(strcat(OutputPath_Epan , '\ScenarioMIP_' , ssp , '\Epan_' , ssp , '_' , GCM) , 'Epan');
    
    clear C r1 Epan
end