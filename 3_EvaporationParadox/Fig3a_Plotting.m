function Fig3a_Plotting(GridEpan_CMIP , Path_Fig3_Output)
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

%% (2) linear regression of AI calculated by [Princeton 1948-2014]
% E_pan m/s to (mm/year)
E_pan_Princeton =  GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_Princeton,1)
    for i_lat = 1 : size(E_pan_Princeton,2)
        if ~isnan(E_pan_Princeton(i_lon , i_lat , 1))
            Y = E_pan_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_E_pan_Princeton_Year(i_lon,i_lat) = b(2);
            p_E_pan_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_E_pan_Princeton_Year(i_lon,i_lat) = nan;
            p_E_pan_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon
k_E_pan_Princeton_Year(isnan(k_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Princeton'],extent,k_E_pan_Princeton_Year');
p_E_pan_Princeton_Year(isnan(p_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Princeton'],extent,p_E_pan_Princeton_Year');
clear k_E_pan_Princeton_Year p_E_pan_Princeton_Year E_pan_Princeton

%% (3) linear regression of AI calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % E_pan m/s to (mm/year)
    E_pan_CMIP =  GridEpan_CMIP(i_ssp+1).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
    for i_lon = 1 : size(E_pan_CMIP,1)
        for i_lat = 1 : size(E_pan_CMIP,2)
            if ~isnan(E_pan_CMIP(i_lon , i_lat , 1))
                Y = E_pan_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_E_pan_CMIP_Year(i_lon,i_lat) = b(2);
                p_E_pan_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_E_pan_CMIP_Year(i_lon,i_lat) = nan;
                p_E_pan_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_E_pan_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_E_pan_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_E_pan_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_E_pan_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_E_pan_CMIP_Year(isnan(k_E_pan_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_' ssp{i_ssp}],extent,k_E_pan_CMIP_Year');
    p_E_pan_CMIP_Year(isnan(p_E_pan_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_' ssp{i_ssp}],extent,p_E_pan_CMIP_Year');
    clear k_E_pan_CMIP_Year p_E_pan_CMIP_Year E_pan_CMIP
end

%% (4) linear regression of AI calculated by CMIP historical experiments
% E_pan m/s to (mm/year)
E_pan_CMIP =  GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_CMIP,1)
    for i_lat = 1 : size(E_pan_CMIP,2)
        if ~isnan(E_pan_CMIP(i_lon , i_lat , 1))
            Y = E_pan_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_E_pan_CMIP_Year(i_lon,i_lat) = b(2);
            p_E_pan_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_E_pan_CMIP_Year(i_lon,i_lat) = nan;
            p_E_pan_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_E_pan_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_E_pan_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_E_pan_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_E_pan_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_E_pan_CMIP_Year(isnan(k_E_pan_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Historical'],extent,k_E_pan_CMIP_Year');
p_E_pan_CMIP_Year(isnan(p_E_pan_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Historical'],extent,p_E_pan_CMIP_Year');
clear k_E_pan_CMIP_Year p_E_pan_CMIP_Year E_pan_CMIP
end
