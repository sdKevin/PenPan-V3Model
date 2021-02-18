function Fig4a_Plotting(GridMet_CMIP , Path_Fig4_Output)
load LandInfo_05deg
%% (1) Adjust map range from 0~360 to -180~180
extent = [-179.75 , 179.75 , -59.75+0.195 , 89.75+0.195];

A = landmask_05deg(1:360 , :); B = landmask_05deg(361:end , :);
landmask_05deg = [B;A]; clear B A
A = lat_05deg(1:360 , :); B = lat_05deg(361:end , :);
lat_05deg = [B;A]; clear B A
A = lon_05deg(1:360 , :); B = lon_05deg(361:end , :) - 360;
lon_05deg = [B;A]; clear B A

for ii = 1 : size(GridMet_CMIP(1).Ensemble_Mean_Met.Sg , 3)
    % Princeton Sg
    A = GridMet_CMIP(1).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(1).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(1).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    % Princeton Li
    A = GridMet_CMIP(1).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(1).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(1).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    % Princeton U10
    A = GridMet_CMIP(1).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(1).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(1).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    % Princeton Ta
    A = GridMet_CMIP(1).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(1).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(1).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    % Princeton Sh
    A = GridMet_CMIP(1).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(1).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(1).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
end
for ii = 1 : size(GridMet_CMIP(2).Ensemble_Mean_Met.Sg , 3)
    % ssp126
    A = GridMet_CMIP(2).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(2).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(2).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(2).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(2).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(2).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(2).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(2).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(2).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(2).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(2).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(2).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(2).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(2).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(2).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
    % ssp245
    A = GridMet_CMIP(3).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(3).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(3).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(3).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(3).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(3).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(3).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(3).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(3).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(3).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(3).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(3).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(3).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(3).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(3).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
    % ssp370
    A = GridMet_CMIP(4).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(4).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(4).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(4).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(4).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(4).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(4).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(4).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(4).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(4).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(4).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(4).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(4).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(4).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(4).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
    % ssp585
    A = GridMet_CMIP(5).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(5).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(5).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(5).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(5).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(5).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(5).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(5).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(5).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(5).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(5).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(5).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(5).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(5).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(5).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
end
for ii = 1 : size(GridMet_CMIP(6).Ensemble_Mean_Met.Sg , 3)
    % Historical
    A = GridMet_CMIP(6).Ensemble_Mean_Met.Sg(1:360 , : , ii);
    B = GridMet_CMIP(6).Ensemble_Mean_Met.Sg(361:end , : , ii);
    GridMet_CMIP(6).Ensemble_Mean_Met.Sg(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(6).Ensemble_Mean_Met.Li(1:360 , : , ii);
    B = GridMet_CMIP(6).Ensemble_Mean_Met.Li(361:end , : , ii);
    GridMet_CMIP(6).Ensemble_Mean_Met.Li(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(6).Ensemble_Mean_Met.U10(1:360 , : , ii);
    B = GridMet_CMIP(6).Ensemble_Mean_Met.U10(361:end , : , ii);
    GridMet_CMIP(6).Ensemble_Mean_Met.U10(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(6).Ensemble_Mean_Met.Ta(1:360 , : , ii);
    B = GridMet_CMIP(6).Ensemble_Mean_Met.Ta(361:end , : , ii);
    GridMet_CMIP(6).Ensemble_Mean_Met.Ta(: , : , ii) = [B;A];
    
    A = GridMet_CMIP(6).Ensemble_Mean_Met.Sh(1:360 , : , ii);
    B = GridMet_CMIP(6).Ensemble_Mean_Met.Sh(361:end , : , ii);
    GridMet_CMIP(6).Ensemble_Mean_Met.Sh(: , : , ii) = [B;A];
end
clear ii A B
clear ii A B

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) linear regression of Met Variables calculated by [Princeton 1948-2014]
% Sg W/m2
Sg_Princeton =  GridMet_CMIP(1).Ensemble_Mean_Met.Sg;
for i_lon = 1 : size(Sg_Princeton,1)
    for i_lat = 1 : size(Sg_Princeton,2)
        if ~isnan(Sg_Princeton(i_lon , i_lat , 1))
            Y = Sg_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Sg_Princeton_Year(i_lon,i_lat) = b(2);
            p_Sg_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Sg_Princeton_Year(i_lon,i_lat) = nan;
            p_Sg_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Sg_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Sg_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Sg_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Sg_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Sg_Princeton_Year(isnan(k_Sg_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Sg_Princeton'],extent,k_Sg_Princeton_Year');
p_Sg_Princeton_Year(isnan(p_Sg_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Sg_Princeton'],extent,p_Sg_Princeton_Year');
clear k_Sg_Princeton_Year p_Sg_Princeton_Year Sg_Princeton

%% (3) linear regression of Met Variables calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % Sg W/m2
    Sg_CMIP =  GridMet_CMIP(i_ssp+1).Ensemble_Mean_Met.Sg;
    for i_lon = 1 : size(Sg_CMIP,1)
        for i_lat = 1 : size(Sg_CMIP,2)
            if ~isnan(Sg_CMIP(i_lon , i_lat , 1))
                Y = Sg_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_Sg_CMIP_Year(i_lon,i_lat) = b(2);
                p_Sg_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_Sg_CMIP_Year(i_lon,i_lat) = nan;
                p_Sg_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_Sg_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_Sg_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_Sg_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_Sg_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_Sg_CMIP_Year(isnan(k_Sg_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'k_Sg_' ssp{i_ssp}],extent,k_Sg_CMIP_Year');
    p_Sg_CMIP_Year(isnan(p_Sg_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'p_Sg_' ssp{i_ssp}],extent,p_Sg_CMIP_Year');
    clear k_Sg_CMIP_Year p_Sg_CMIP_Year Sg_CMIP
end

%% (4) linear regression of Met Variables calculated by CMIP historical experiments
% Sg W/m2
Sg_CMIP =  GridMet_CMIP(6).Ensemble_Mean_Met.Sg;
for i_lon = 1 : size(Sg_CMIP,1)
    for i_lat = 1 : size(Sg_CMIP,2)
        if ~isnan(Sg_CMIP(i_lon , i_lat , 1))
            Y = Sg_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Sg_CMIP_Year(i_lon,i_lat) = b(2);
            p_Sg_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Sg_CMIP_Year(i_lon,i_lat) = nan;
            p_Sg_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Sg_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Sg_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Sg_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Sg_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Sg_CMIP_Year(isnan(k_Sg_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Sg_Historical'],extent,k_Sg_CMIP_Year');
p_Sg_CMIP_Year(isnan(p_Sg_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Sg_Historical'],extent,p_Sg_CMIP_Year');
clear k_Sg_CMIP_Year p_Sg_CMIP_Year Sg_CMIP


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Li %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) linear regression of Met Variables calculated by [Princeton 1948-2014]
% Li W/m2
Li_Princeton =  GridMet_CMIP(1).Ensemble_Mean_Met.Li;
for i_lon = 1 : size(Li_Princeton,1)
    for i_lat = 1 : size(Li_Princeton,2)
        if ~isnan(Li_Princeton(i_lon , i_lat , 1))
            Y = Li_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Li_Princeton_Year(i_lon,i_lat) = b(2);
            p_Li_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Li_Princeton_Year(i_lon,i_lat) = nan;
            p_Li_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Li_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Li_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Li_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Li_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Li_Princeton_Year(isnan(k_Li_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Li_Princeton'],extent,k_Li_Princeton_Year');
p_Li_Princeton_Year(isnan(p_Li_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Li_Princeton'],extent,p_Li_Princeton_Year');
clear k_Li_Princeton_Year p_Li_Princeton_Year Li_Princeton

%% (3) linear regression of Met Variables calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % Li W/m2
    Li_CMIP =  GridMet_CMIP(i_ssp+1).Ensemble_Mean_Met.Li;
    for i_lon = 1 : size(Li_CMIP,1)
        for i_lat = 1 : size(Li_CMIP,2)
            if ~isnan(Li_CMIP(i_lon , i_lat , 1))
                Y = Li_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_Li_CMIP_Year(i_lon,i_lat) = b(2);
                p_Li_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_Li_CMIP_Year(i_lon,i_lat) = nan;
                p_Li_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_Li_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_Li_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_Li_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_Li_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_Li_CMIP_Year(isnan(k_Li_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'k_Li_' ssp{i_ssp}],extent,k_Li_CMIP_Year');
    p_Li_CMIP_Year(isnan(p_Li_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'p_Li_' ssp{i_ssp}],extent,p_Li_CMIP_Year');
    clear k_Li_CMIP_Year p_Li_CMIP_Year Li_CMIP
end

%% (4) linear regression of Met Variables calculated by CMIP historical experiments
% Li W/m2
Li_CMIP =  GridMet_CMIP(6).Ensemble_Mean_Met.Li;
for i_lon = 1 : size(Li_CMIP,1)
    for i_lat = 1 : size(Li_CMIP,2)
        if ~isnan(Li_CMIP(i_lon , i_lat , 1))
            Y = Li_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Li_CMIP_Year(i_lon,i_lat) = b(2);
            p_Li_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Li_CMIP_Year(i_lon,i_lat) = nan;
            p_Li_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Li_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Li_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Li_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Li_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Li_CMIP_Year(isnan(k_Li_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Li_Historical'],extent,k_Li_CMIP_Year');
p_Li_CMIP_Year(isnan(p_Li_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Li_Historical'],extent,p_Li_CMIP_Year');
clear k_Li_CMIP_Year p_Li_CMIP_Year Li_CMIP


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% U10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) linear regression of Met Variables calculated by [Princeton 1948-2014]
% U10 m/s
U10_Princeton =  GridMet_CMIP(1).Ensemble_Mean_Met.U10;
for i_lon = 1 : size(U10_Princeton,1)
    for i_lat = 1 : size(U10_Princeton,2)
        if ~isnan(U10_Princeton(i_lon , i_lat , 1))
            Y = U10_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_U10_Princeton_Year(i_lon,i_lat) = b(2);
            p_U10_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_U10_Princeton_Year(i_lon,i_lat) = nan;
            p_U10_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_U10_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_U10_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_U10_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_U10_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_U10_Princeton_Year(isnan(k_U10_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_U10_Princeton'],extent,k_U10_Princeton_Year');
p_U10_Princeton_Year(isnan(p_U10_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_U10_Princeton'],extent,p_U10_Princeton_Year');
clear k_U10_Princeton_Year p_U10_Princeton_Year U10_Princeton

%% (3) linear regression of Met Variables calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % U10 m/s
    U10_CMIP =  GridMet_CMIP(i_ssp+1).Ensemble_Mean_Met.U10;
    for i_lon = 1 : size(U10_CMIP,1)
        for i_lat = 1 : size(U10_CMIP,2)
            if ~isnan(U10_CMIP(i_lon , i_lat , 1))
                Y = U10_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_U10_CMIP_Year(i_lon,i_lat) = b(2);
                p_U10_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_U10_CMIP_Year(i_lon,i_lat) = nan;
                p_U10_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_U10_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_U10_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_U10_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_U10_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_U10_CMIP_Year(isnan(k_U10_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'k_U10_' ssp{i_ssp}],extent,k_U10_CMIP_Year');
    p_U10_CMIP_Year(isnan(p_U10_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'p_U10_' ssp{i_ssp}],extent,p_U10_CMIP_Year');
    clear k_U10_CMIP_Year p_U10_CMIP_Year U10_CMIP
end

%% (4) linear regression of Met Variables calculated by CMIP historical experiments
% U10 m/s
U10_CMIP =  GridMet_CMIP(6).Ensemble_Mean_Met.U10;
for i_lon = 1 : size(U10_CMIP,1)
    for i_lat = 1 : size(U10_CMIP,2)
        if ~isnan(U10_CMIP(i_lon , i_lat , 1))
            Y = U10_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_U10_CMIP_Year(i_lon,i_lat) = b(2);
            p_U10_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_U10_CMIP_Year(i_lon,i_lat) = nan;
            p_U10_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_U10_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_U10_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_U10_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_U10_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_U10_CMIP_Year(isnan(k_U10_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_U10_Historical'],extent,k_U10_CMIP_Year');
p_U10_CMIP_Year(isnan(p_U10_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_U10_Historical'],extent,p_U10_CMIP_Year');
clear k_U10_CMIP_Year p_U10_CMIP_Year U10_CMIP


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) linear regression of Met Variables calculated by [Princeton 1948-2014]
% Ta k
Ta_Princeton =  GridMet_CMIP(1).Ensemble_Mean_Met.Ta;
for i_lon = 1 : size(Ta_Princeton,1)
    for i_lat = 1 : size(Ta_Princeton,2)
        if ~isnan(Ta_Princeton(i_lon , i_lat , 1))
            Y = Ta_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Ta_Princeton_Year(i_lon,i_lat) = b(2);
            p_Ta_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Ta_Princeton_Year(i_lon,i_lat) = nan;
            p_Ta_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Ta_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Ta_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Ta_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Ta_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Ta_Princeton_Year(isnan(k_Ta_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Ta_Princeton'],extent,k_Ta_Princeton_Year');
p_Ta_Princeton_Year(isnan(p_Ta_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Ta_Princeton'],extent,p_Ta_Princeton_Year');
clear k_Ta_Princeton_Year p_Ta_Princeton_Year Ta_Princeton

%% (3) linear regression of Met Variables calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % Ta K
    Ta_CMIP =  GridMet_CMIP(i_ssp+1).Ensemble_Mean_Met.Ta;
    for i_lon = 1 : size(Ta_CMIP,1)
        for i_lat = 1 : size(Ta_CMIP,2)
            if ~isnan(Ta_CMIP(i_lon , i_lat , 1))
                Y = Ta_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_Ta_CMIP_Year(i_lon,i_lat) = b(2);
                p_Ta_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_Ta_CMIP_Year(i_lon,i_lat) = nan;
                p_Ta_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_Ta_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_Ta_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_Ta_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_Ta_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_Ta_CMIP_Year(isnan(k_Ta_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'k_Ta_' ssp{i_ssp}],extent,k_Ta_CMIP_Year');
    p_Ta_CMIP_Year(isnan(p_Ta_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'p_Ta_' ssp{i_ssp}],extent,p_Ta_CMIP_Year');
    clear k_Ta_CMIP_Year p_Ta_CMIP_Year Ta_CMIP
end

%% (4) linear regression of Met Variables calculated by CMIP historical experiments
% Ta K
Ta_CMIP =  GridMet_CMIP(6).Ensemble_Mean_Met.Ta;
for i_lon = 1 : size(Ta_CMIP,1)
    for i_lat = 1 : size(Ta_CMIP,2)
        if ~isnan(Ta_CMIP(i_lon , i_lat , 1))
            Y = Ta_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Ta_CMIP_Year(i_lon,i_lat) = b(2);
            p_Ta_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Ta_CMIP_Year(i_lon,i_lat) = nan;
            p_Ta_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Ta_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Ta_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Ta_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Ta_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Ta_CMIP_Year(isnan(k_Ta_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Ta_Historical'],extent,k_Ta_CMIP_Year');
p_Ta_CMIP_Year(isnan(p_Ta_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Ta_Historical'],extent,p_Ta_CMIP_Year');
clear k_Ta_CMIP_Year p_Ta_CMIP_Year Ta_CMIP


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) linear regression of Met Variables calculated by [Princeton 1948-2014]
% Sh
Sh_Princeton =  GridMet_CMIP(1).Ensemble_Mean_Met.Sh;
for i_lon = 1 : size(Sh_Princeton,1)
    for i_lat = 1 : size(Sh_Princeton,2)
        if ~isnan(Sh_Princeton(i_lon , i_lat , 1))
            Y = Sh_Princeton(i_lon , i_lat , :);
            Y = Y(:);
            x = [1948:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Sh_Princeton_Year(i_lon,i_lat) = b(2);
            p_Sh_Princeton_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Sh_Princeton_Year(i_lon,i_lat) = nan;
            p_Sh_Princeton_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Sh_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Sh_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Sh_Princeton_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Sh_Princeton_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Sh_Princeton_Year(isnan(k_Sh_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Sh_Princeton'],extent,k_Sh_Princeton_Year');
p_Sh_Princeton_Year(isnan(p_Sh_Princeton_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Sh_Princeton'],extent,p_Sh_Princeton_Year');
clear k_Sh_Princeton_Year p_Sh_Princeton_Year Sh_Princeton

%% (3) linear regression of Met Variables calculated by CMIP scenarios
ssp = {'ssp126','ssp245','ssp370','ssp585'};
for i_ssp = 1 : length(ssp)
    % Sh
    Sh_CMIP =  GridMet_CMIP(i_ssp+1).Ensemble_Mean_Met.Sh;
    for i_lon = 1 : size(Sh_CMIP,1)
        for i_lat = 1 : size(Sh_CMIP,2)
            if ~isnan(Sh_CMIP(i_lon , i_lat , 1))
                Y = Sh_CMIP(i_lon , i_lat , :);
                Y = Y(:);
                x = [2015:2100];
                X = [ones(length(Y),1) , x'];
                [b,bint,r,rint,stats] = regress(Y,X);
                k_Sh_CMIP_Year(i_lon,i_lat) = b(2);
                p_Sh_CMIP_Year(i_lon,i_lat) = stats(3);
                clear x Y X b bint r rint stats
            else
                k_Sh_CMIP_Year(i_lon,i_lat) = nan;
                p_Sh_CMIP_Year(i_lon,i_lat) = nan;
            end
        end
    end
    clear i_lat i_lon
    
    % interpolate the seam
    k_Sh_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),k_Sh_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    p_Sh_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
        lon_05deg([1:358,362:end],:),p_Sh_CMIP_Year([1:358,362:end],:),...
        lat_05deg,lon_05deg).*landmask_05deg;
    
    k_Sh_CMIP_Year(isnan(k_Sh_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'k_Sh_' ssp{i_ssp}],extent,k_Sh_CMIP_Year');
    p_Sh_CMIP_Year(isnan(p_Sh_CMIP_Year)) = -9999;
    SaveData2GeoTIFF([Path_Fig4_Output 'p_Sh_' ssp{i_ssp}],extent,p_Sh_CMIP_Year');
    clear k_Sh_CMIP_Year p_Sh_CMIP_Year Sh_CMIP
end

%% (4) linear regression of Met Variables calculated by CMIP historical experiments
% Sh
Sh_CMIP =  GridMet_CMIP(6).Ensemble_Mean_Met.Sh;
for i_lon = 1 : size(Sh_CMIP,1)
    for i_lat = 1 : size(Sh_CMIP,2)
        if ~isnan(Sh_CMIP(i_lon , i_lat , 1))
            Y = Sh_CMIP(i_lon , i_lat , :);
            Y = Y(:);
            x = [1850:2014];
            X = [ones(length(Y),1) , x'];
            [b,bint,r,rint,stats] = regress(Y,X);
            k_Sh_CMIP_Year(i_lon,i_lat) = b(2);
            p_Sh_CMIP_Year(i_lon,i_lat) = stats(3);
            clear x Y X b bint r rint stats
        else
            k_Sh_CMIP_Year(i_lon,i_lat) = nan;
            p_Sh_CMIP_Year(i_lon,i_lat) = nan;
        end
    end
end
clear i_lat i_lon

% interpolate the seam
k_Sh_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),k_Sh_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;
p_Sh_CMIP_Year = interp2(lat_05deg([1:358,362:end],:),...
    lon_05deg([1:358,362:end],:),p_Sh_CMIP_Year([1:358,362:end],:),...
    lat_05deg,lon_05deg).*landmask_05deg;

k_Sh_CMIP_Year(isnan(k_Sh_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'k_Sh_Historical'],extent,k_Sh_CMIP_Year');
p_Sh_CMIP_Year(isnan(p_Sh_CMIP_Year)) = -9999;
SaveData2GeoTIFF([Path_Fig4_Output 'p_Sh_Historical'],extent,p_Sh_CMIP_Year');
clear k_Sh_CMIP_Year p_Sh_CMIP_Year Sh_CMIP

end
