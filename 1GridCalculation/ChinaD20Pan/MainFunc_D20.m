clc; clear all; close all;

%% (1) Setting the input/output paths
% CMIP6 Ensemble Meteorological Data
InputPath_CMIP6_Ensemble = 'D:\CMIP6\ProcessData\Ensemble_Met';
% Princeton-GMFD Data
InputPath_Princeton = 'D:\CMIP6\ProcessData\Princeton\monthly';
% Save Meteorological variables
OutputPath_MetVar = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_Met';
% Save Pan evaporation (Epan) Data
OutputPath_ETpan = 'E:\PenPanV3\VariableStorage\MonthlyVar\Var_D20Epan';

%% (2) Historical Experiment
% Name of Global Climate Model
GCM_Ensemble = {'ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-CanOE',...
    'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg',...
    'FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8',...
    'INM-CM5-0','IPSL-CM6A-LR','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
    'MRI-ESM2-0','NorESM2-MM','UKESM1-0-LL'};
for i_GCM = 1 : length(GCM_Ensemble)
    %% (2.1) load Data
    GCM = GCM_Ensemble{i_GCM}
    load(strcat(InputPath_CMIP6_Ensemble , '\Historical\' , GCM , '.mat'));
    
    %% (2.2) Interpolating Forcing Data to Uniform Resolution i.e., 0.5deg
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
        R1.sfcWind(:,:,ii) = interp2(r1.lat , r1.lon , r1.sfcWind(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
        R1.tas(:,:,ii) = interp2(r1.lat , r1.lon , r1.tas(:,:,ii) ,...
            lat_05deg , lon_05deg) .* landmask_05deg;
    end
    clear r1 ii elevation_05deg
    r1 = R1; clear R1
    
    %% (2.3) Bias Correction : Using Princeton-GMFD Data to correct CMIP6 Historical Data
    load([InputPath_Princeton , '\huss.mat']); huss(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.huss = nanmean(huss,3) ./ nanmean(r1.huss(:,:,1177:end),3); C.huss(C.huss>10)=10;clear huss;
    load([InputPath_Princeton , '\ps.mat']); ps(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.ps = nanmean(ps,3) ./ nanmean(r1.ps(:,:,1177:end),3); C.ps(C.ps>10)=10;clear ps;
    load([InputPath_Princeton , '\rlds.mat']); rlds(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.rlds = nanmean(rlds,3) ./ nanmean(r1.rlds(:,:,1177:end),3); C.rlds(C.rlds>10)=10;clear rlds;
    load([InputPath_Princeton , '\rsds.mat']); rsds(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.rsds = nanmean(rsds,3) ./ nanmean(r1.rsds(:,:,1177:end),3); C.rsds(C.rsds>10)=10;clear rsds;
    load([InputPath_Princeton , '\sfcWind.mat']); sfcWind(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.sfcWind = nanmean(sfcWind,3) ./ nanmean(r1.sfcWind(:,:,1177:end),3); C.sfcWind(C.sfcWind>10)=10;clear sfcWind;
    load([InputPath_Princeton , '\tas.mat']); tas(:,:,805:828)=[]; %1948-2016 to 1948-2014
    C.tas = nanmean(tas,3) - nanmean(r1.tas(:,:,1177:end),3);clear tas;
    for ii = 1:size(r1.huss,3)
        r1.huss(:,:,ii) = r1.huss(:,:,ii) .* C.huss;
        r1.ps(:,:,ii) = r1.ps(:,:,ii) .* C.ps;
        r1.rlds(:,:,ii) = r1.rlds(:,:,ii) .* C.rlds;
        r1.rsds(:,:,ii) = r1.rsds(:,:,ii) .* C.rsds;
        r1.sfcWind(:,:,ii) = r1.sfcWind(:,:,ii) .* C.sfcWind;
        r1.tas(:,:,ii) = r1.tas(:,:,ii) + C.tas;
    end
    clear ii
    
    %% (2.4) Calculation pan evaporation (Epan) (m/s)
    %% 2 计算蒸发皿的参数 (D20)
    % L is the diameter of Class A
    pan_pars.D = 0.2; % [m]
    pan_pars.L = pan_pars.D; % [m]
    % he is the height of rim
    pan_pars.he = 0.08; % [m]
    % hw is the height of water level
    pan_pars.hw = 0.02; % [m]
    % Beta is the ratio of heat to mass transfer coefficients of the pan
    pan_pars.Beta = 2 + 0.2*0.1./(0.25*pi*pan_pars.D^2) + pi*pan_pars.D*0.1./(0.25*pi*pan_pars.D^2) + 0.2*0.08./(0.25*pi*pan_pars.D^2) + pi*pan_pars.D*0.08./(0.25*pi*pan_pars.D^2); % Class A 0.97 外壁面积 0.31 断面面积 1.15 水面面积
    % C is the correction factor to account for the shading effect of the bird guard
    pan_pars.C = 1; %Class A C=1.07; D20 C=1; 601B C=1;
    % e_gnd is the emissivity of ground
    pan_pars.e_gnd = 0.90; % LWH2013 这个数可以保持不变
    % e_wall is the emissivity of water
    pan_pars.e_w = 0.89; % LWH2013 这个数可以保持不变
    % e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
    pan_pars.e_wall = 0.82; %D20可以使用类似的数据（同一个参考文献）因为外表同为ZnO薄膜，但是601B需要重新确定数值。
    % is the refractive index
    pan_pars.N = 1.33; % Class A: 1.33
    % is the extinction coeeficient
    pan_pars.K = 0; % Class A: 0
    % The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
    pan_pars.alpha_0_wall = 0.36;
    % alpha_gnd
    pan_pars.alpha_gnd = 0.2;
    
    Epan = PenPan_V3_D20(pan_pars_n , pan_pars_q, pan_pars , lat_05deg , elevation_05deg ,...
        r1.rsds , r1.rlds , r1.sfcWind , r1.tas , r1.huss , r1.ps ,...
        [OutputPath_MetVar , '\Historical'] , strcat('Historical_',GCM));
    
    
    
    
    
    
    
    
    
    PM_RC = Penman_Mothith(r1.rsds , r1.rlds , r1.sfcWind , r1.tas , r1.huss , r1.ps);
    [PM_RC_CO2_Yang , PM_RC_CO2_Jarvis_H , PM_RC_CO2_Jarvis_L] =...
        Penman_Mothith_CO2(r1.rsds , r1.rlds , r1.sfcWind , r1.tas , r1.huss , r1.ps , CO2, r1.pr ,...
        [OutputPath_Attribution , '\Historical'] , [OutputPath_MetVar , '\Historical'] , strcat('Historical_',GCM));
    
    %% (2.5) Save the result
    % strcat(OutputPath_ETrc , '\Historical\ETrc_',GCM) is for saving 'PM_RC','PM_RC_CO2_Yang','PM_RC_CO2_Jarvis_H','PM_RC_CO2_Jarvis_L' in
    % ETrc_GCM.mat, which is the reference crop evapotranspiration calculated by different methods
    save(strcat(OutputPath_ETrc , '\Historical\ETrc_Historical_',GCM),...
        'PM_RC','PM_RC_CO2_Yang','PM_RC_CO2_Jarvis_H','PM_RC_CO2_Jarvis_L');
    
    clear C CO2 r1 PM_RC PM_RC_CO2_Yang PM_RC_CO2_Jarvis_H PM_RC_CO2_Jarvis_L
end












%% 1 DataPrep
[ dem , Tmax_Month , Tmin_Month , U2_Month , SH_Month ,...
    RH_Month , Mask , pars_Efv_D20_n , pars_Efv_D20_q , longd_Mon , latd_Mon , as, bs] = DataPrep();
% 'Latitude(°),','Longitude(°),','Elevation,'（m）
% 'RHU,'(%)【没有均一化】,'SSD,'(1hr),'MaxT,'(℃),'MinT,'(℃),'Win'(m/s),as,bs
%% 3 计算蒸发量
for Year = 1960:2014
    for Month = 1:12
        XuHao = (Year-1960).*12 + Month;
        if Month<10
            Date = [num2str(Year) '0' num2str(Month)];
            Date = str2num(Date);
        else
            Date = [num2str(Year) num2str(Month)];
            Date = str2num(Date);
        end
        [Eva_Cal(:,:,XuHao),Rn_lamda_rou(:,:,XuHao),E_pan_A(:,:,XuHao),E_pan_R(:,:,XuHao),...
            E_pan_R_w(:,:,XuHao), E_pan_R_wall(:,:,XuHao),E_pan_R_rim(:,:,XuHao) ,...
            E_pan_R_rim_pie(:,:,XuHao),E_pan_R_bot(:,:,XuHao)]...
            = PenPan_V2_Steady_D20(pars_Efv_D20_n(:,:), pars_Efv_D20_q(:,:) , pan_pars, latd_Mon(:,:) , dem(:,:) ,...
            RH_Month(:,:,XuHao) , Tmax_Month(:,:,XuHao) , Tmin_Month(:,:,XuHao) , U2_Month(:,:,XuHao) ,...
            SH_Month(:,:,XuHao) , as(:,:,Month) , bs(:,:,Month) , Date);
        
        Eva_Cal(:,:,XuHao) = Eva_Cal(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        Rn_lamda_rou(:,:,XuHao) = Rn_lamda_rou(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_A(:,:,XuHao) = E_pan_A(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R(:,:,XuHao) = E_pan_R(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R_w(:,:,XuHao) = E_pan_R_w(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R_wall(:,:,XuHao) = E_pan_R_wall(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R_rim(:,:,XuHao) = E_pan_R_rim(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R_rim_pie(:,:,XuHao) = E_pan_R_rim_pie(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        E_pan_R_bot(:,:,XuHao) = E_pan_R_bot(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s 转换成 mm
        %运行所有月份需要2天时间
    end
    Year
end
Eva_Cal(Eva_Cal<0)=0;
Rn_lamda_rou(Rn_lamda_rou<0)=0;
E_pan_A(E_pan_A<0)=0;
E_pan_R(E_pan_R<0)=0;
E_pan_R_w(E_pan_R_w<0)=0;
E_pan_R_wall(E_pan_R_wall<0)=0;
E_pan_R_rim(E_pan_R_rim<0)=0;
E_pan_R_rim_pie(E_pan_R_rim_pie<0)=0;
E_pan_R_bot(E_pan_R_bot<0)=0;
save Eva_Cal Eva_Cal
save E_pan_A E_pan_A
save E_pan_R E_pan_R
save E_pan_R_w E_pan_R_w
save E_pan_R_wall E_pan_R_wall
save E_pan_R_rim E_pan_R_rim
save E_pan_R_rim_pie E_pan_R_rim_pie
save E_pan_R_bot E_pan_R_bot
%% 4 画图
Plotting( Eva_Cal )
%% 5 保存为nc数据
% 保存 Eva_Cal
Value_lon = [70.025:0.05:140.025];
Value_lat = [15.025:0.05:55.025];
Out_Path = ['D20_Epan.nc'];
Object='D20_Epan';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Eva_Cal,Object,days_since,days_end);
% ncdisp(['D20_Epan.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['D20_Epan.nc'], Object);% read variables in nc file
% 保存 Epan_wb
Epan_wb = E_pan_A + E_pan_R_w;
Out_Path = ['D20_Epan_wb.nc'];
Object='D20_Epan_wb';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Epan_wb,Object,days_since,days_end);
% ncdisp(['D20_Epan_wb.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['D20_Epan_wb.nc'], Object);% read variables in nc file
% 保存 Epan_pw
Epan_pw = E_pan_R - E_pan_R_w;
Out_Path = ['D20_Epan_pw.nc'];
Object='D20_Epan_pw';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Epan_pw,Object,days_since,days_end);
% ncdisp(['D20_Epan_pw.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['D20_Epan_pw.nc'], Object);% read variables in nc file