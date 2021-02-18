function Fig5a_Plotting(GridEpan_CMIP , Path_Fig5_Output)
load LandInfo_05deg
%% (1) Adjust map range from 0~360 to -180~180
extent = [-179.75 , 179.75 , -59.75+0.195 , 89.75+0.195];

A = landmask_05deg(1:360 , :); B = landmask_05deg(361:end , :);
landmask_05deg = [B;A]; clear B A
A = lat_05deg(1:360 , :); B = lat_05deg(361:end , :);
lat_05deg = [B;A]; clear B A
A = lon_05deg(1:360 , :); B = lon_05deg(361:end , :) - 360;
lon_05deg = [B;A]; clear B A

for ii = 1 : size(GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan , 3)
    % Princeton E_pan
    A = GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
end
for ii = 1 : size(GridEpan_CMIP(2).Ensemble_Mean_Epan.E_pan , 3)
    % ssp126 E_pan
    A = GridEpan_CMIP(2).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(2).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(2).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
    % ssp245 E_pan
    A = GridEpan_CMIP(3).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(3).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(3).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
    % ssp370 E_pan
    A = GridEpan_CMIP(4).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(4).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(4).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
    % ssp585 E_pan
    A = GridEpan_CMIP(5).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(5).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(5).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
end
for ii = 1 : size(GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan , 3)
    % Historical E_pan
    A = GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
end
clear ii A B
clear ii A B

%% (2) Epan calculated by [Princeton 1948-2014]
% E_pan m/s to (mm/year)
E_pan_Princeton =  GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
E_pan_Mean_Princeton = nanmean(E_pan_Princeton,3);

% interpolate the seam
E_pan_Mean_Princeton = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),E_pan_Mean_Princeton([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

E_pan_Mean_Princeton(isnan(E_pan_Mean_Princeton)) = -9999;
SaveData2GeoTIFF([Path_Fig5_Output 'E_pan_Mean_Princeton'],extent,E_pan_Mean_Princeton');
clear E_pan_Mean_Princeton E_pan_Princeton

%% (3) Epan calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % E_pan m/s to (mm/year)
    E_pan_CMIP =  GridEpan_CMIP(i_ssp+1).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
    E_pan_Mean_CMIP = nanmean(E_pan_CMIP,3);
    
    % interpolate the seam
    E_pan_Mean_CMIP = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),E_pan_Mean_CMIP([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    E_pan_Mean_CMIP(isnan(E_pan_Mean_CMIP)) = -9999;
    SaveData2GeoTIFF([Path_Fig5_Output 'E_pan_Mean_' ssp{i_ssp}],extent,E_pan_Mean_CMIP');
    clear E_pan_Mean_CMIP E_pan_CMIP
end

%% (4) linear regression of Epan calculated by CMIP historical experiments
% E_pan m/s to (mm/year)
E_pan_CMIP =  GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
E_pan_Mean_CMIP = nanmean(E_pan_CMIP,3);

% interpolate the seam
E_pan_Mean_CMIP = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),E_pan_Mean_CMIP([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

E_pan_Mean_CMIP(isnan(E_pan_Mean_CMIP)) = -9999;
SaveData2GeoTIFF([Path_Fig5_Output 'E_pan_Mean_Historical'],extent,E_pan_Mean_CMIP');
clear E_pan_Mean_CMIP E_pan_CMIP
end
