function Fig3b_Plotting(GridEpan_CMIP , Path_Fig3_Output)
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
for ii = 1 : size(GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan , 3)
    % Historical E_pan
    A = GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(1:360 , : , ii);
    B = GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(361:end , : , ii);
    GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(: , : , ii) = [B;A];
end
clear ii A B
clear ii A B

%% (2) linear regression of Epan calculated by [Princeton 1948-1992]
% E_pan m/s to (mm/year)
E_pan_Princeton =  GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan(:,:,1:45) .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_Princeton,1)
    for i_lat = 1 : size(E_pan_Princeton,2)
        if ~isnan(E_pan_Princeton(i_lon , i_lat , 1))
            Y = E_pan_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:1992];
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

% interpolate the seam
k_E_pan_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_E_pan_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_E_pan_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_E_pan_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_E_pan_Princeton_Year(isnan(k_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Princeton_1948_1992'],extent,k_E_pan_Princeton_Year');
p_E_pan_Princeton_Year(isnan(p_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Princeton_1948_1992'],extent,p_E_pan_Princeton_Year');
clear k_E_pan_Princeton_Year p_E_pan_Princeton_Year E_pan_Princeton

%% (3) linear regression of Epan calculated by [Princeton 1993-2014]
% E_pan m/s to (mm/year)
E_pan_Princeton =  GridEpan_CMIP(1).Ensemble_Mean_Epan.E_pan(:,:,46:end) .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_Princeton,1)
    for i_lat = 1 : size(E_pan_Princeton,2)
        if ~isnan(E_pan_Princeton(i_lon , i_lat , 1))
            Y = E_pan_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1993:2014];
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

% interpolate the seam
k_E_pan_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_E_pan_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_E_pan_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_E_pan_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_E_pan_Princeton_Year(isnan(k_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Princeton_1993_2014'],extent,k_E_pan_Princeton_Year');
p_E_pan_Princeton_Year(isnan(p_E_pan_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Princeton_1993_2014'],extent,p_E_pan_Princeton_Year');
clear k_E_pan_Princeton_Year p_E_pan_Princeton_Year E_pan_Princeton

%% (4) linear regression of Epan calculated by CMIP historical experiments (1948-1992)
% E_pan m/s to (mm/year)
E_pan_CMIP =  GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(:,:,99:143) .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_CMIP,1)
    for i_lat = 1 : size(E_pan_CMIP,2)
        if ~isnan(E_pan_CMIP(i_lon , i_lat , 1))
            Y = E_pan_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:1992];
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
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Historical_1948_1992'],extent,k_E_pan_CMIP_Year');
p_E_pan_CMIP_Year(isnan(p_E_pan_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Historical_1948_1992'],extent,p_E_pan_CMIP_Year');
clear k_E_pan_CMIP_Year p_E_pan_CMIP_Year E_pan_CMIP

%% (5) linear regression of Epan calculated by CMIP historical experiments (1993-2014)
% E_pan m/s to (mm/year)
E_pan_CMIP =  GridEpan_CMIP(6).Ensemble_Mean_Epan.E_pan(:,:,144:end) .* 365.*24.*3600.*1000;
for i_lon = 1 : size(E_pan_CMIP,1)
    for i_lat = 1 : size(E_pan_CMIP,2)
        if ~isnan(E_pan_CMIP(i_lon , i_lat , 1))
            Y = E_pan_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1993:2014];
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
SaveData2GeoTIFF([Path_Fig3_Output 'k_E_pan_Historical_1993_2014'],extent,k_E_pan_CMIP_Year');
p_E_pan_CMIP_Year(isnan(p_E_pan_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig3_Output 'p_E_pan_Historical_1993_2014'],extent,p_E_pan_CMIP_Year');
clear k_E_pan_CMIP_Year p_E_pan_CMIP_Year E_pan_CMIP

end
