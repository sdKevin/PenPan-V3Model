clc; clear all; close all;

%% Setting the input/output paths
% Princeton Global Meteorological Forcing (Princeton-GMFD) Data
InputPath_Princeton = 'D:\CMIP6\ProcessData\Princeton\monthly';
% Save Pan evaporation (Epan) Data
OutputPath_Epan = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_601BEpan';

%% Princeton-GMFD Data Preparation
% Load Global 0.5 Degree Coordinate Data from Princeton-GMFD
load LandInfo_05deg
% Load Princeton-GMFD
load([InputPath_Princeton , '\huss.mat']); huss = huss .* landmask_05deg; huss(:,:,805:828) = []; %1948-2016 to 1948-2014
load([InputPath_Princeton , '\ps.mat']); ps = ps .* landmask_05deg; ps(:,:,805:828) = []; %1948-2016 to 1948-2014
load([InputPath_Princeton , '\rlds.mat']); rlds = rlds .* landmask_05deg; rlds(:,:,805:828) = []; %1948-2016 to 1948-2014
load([InputPath_Princeton , '\rsds.mat']); rsds = rsds .* landmask_05deg; rsds(:,:,805:828) = []; %1948-2016 to 1948-2014
load([InputPath_Princeton , '\sfcWind.mat']); sfcWind = sfcWind .* landmask_05deg; sfcWind(:,:,805:828) = []; %1948-2016 to 1948-2014
load([InputPath_Princeton , '\tas.mat']); tas = tas .* landmask_05deg; tas(:,:,805:828) = []; %1948-2016 to 1948-2014
% Load a CMIP6 rsdt (since CMIP estimated rsdt is closed to actual value)
load('E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met\Historical\Met_Var_Historical_BCC-CSM2-MR.mat');
rsdt = Met_Var.Ra(:,:,1177:end); clear Met_Var
%% Calculate pan evaporation (Epan) (m/s)
%% Pan Parameters
% L is the diameter of Class A
pan_pars.D = 0.618; % [m]
pan_pars.L = pan_pars.D + 0.4; % [m] Diameter of outside ring, since which determines the height of boundary layer
% he is the height of rim
pan_pars.he = 0.075; % [m]
% hw is the height of water level
pan_pars.hw = 0.612; % [m]
% Beta is the ratio of heat to mass transfer coefficients of the pan
pan_pars.Beta = 1 + 2*pi*pan_pars.D*pan_pars.he./(0.25*pi*pan_pars.D^2);
% C is the correction factor to account for the shading effect of the bird guard
pan_pars.C = 1; %Class A C=1.07; D20 C=1; 601B C=1;
% e_gnd is the emissivity of ground
pan_pars.e_gnd = 0.90;
% e_wall is the emissivity of water
pan_pars.e_w = 0.89;
% e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
pan_pars.e_wall = 0.85;
% N is the refractive index
pan_pars.N = 1.33; % Class A: 1.33
% K is the extinction coeeficient
pan_pars.K = 0; % Class A: 0
% The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
pan_pars.alpha_0_wall = 0.80;
% alpha_gnd
pan_pars.alpha_gnd = 0.2;

Epan = PenPan_V3_601B(pan_pars , lat_05deg , elevation_05deg ,...
    rsds , rsdt , rlds , sfcWind , tas , huss , ps); % Epan(m/s)

%% Save the result
save(strcat(OutputPath_Epan , '\Princeton\Epan_Princeton') , 'Epan');

clear Epan