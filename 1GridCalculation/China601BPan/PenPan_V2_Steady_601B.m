function [Result,Rn_lamda_rou,E_pan_A,E_pan_R,E_pan_R_w,E_pan_R_rim,E_pan_R_rim_pie]...
    = PenPan_V2_Steady_601B(pan_pars_n , pan_pars_q, pan_pars, latd , Elevation , rhum , Tmax , Tmin , U2 ,...
    sunshinehour , as , bs , Date)
% Lim Wee Ho, 2016, PenPan_V2_Steady
% latd °; Elevation m; rhum %; Tmax ℃; Tmin ℃; U10 m/s; days ; Date 199301-201612
%% 1 单位转变
latd = pi/180 .* latd; %°转rad latd=-17.95*pi/180
Tmax = Tmax + 273.15; %from ℃ to K
Tmin = Tmin + 273.15; %from ℃ to K
%% 2.1 Constants
% Mw [kg/mol] is the molecular mass of water (=0.018 kg/mol).
Mw = 0.018; % Appendix B
% Ma [kg/mol] is the molecular mass of air .
Ma = 0.029;
% R [ J/(mol*K) ] is the ideal gas constant
R = 8.314;
% rou_w [kg/m^3] is the density of liquid water
rou_w = 1000;
% Po [Pa] is the atmospheric pressure at the mean sea level
Po = 101.3 .*1000;
% Sigma [W/(m^2K^4)] is the Stefan-Boltzmann constant
Sigma = 5.67./100000000;
%% 2.2
% 风速函数会用到
% Ta is the air temprature[K]
Ta = (Tmax + Tmin) ./ 2;
% lamda [J/kg] is the latent heat of vaporisation of liquid water
lamda = 2.501 .* 1000000 - 2370 .* ( Ta - 273.15 ); % Appendix B
% Gama [Pa/K] psychrometric constant  [Pa K1]
Gama = 67 - 0.0072 .* Elevation; % Appendix B
% s [Pa/K] is the slope of the saturation vapour pressure versus temperature curve at air temperature (Ta, [K])
s = 611 .* lamda .* Mw .* ( exp( 17.27.*(Ta-273.15)./(Ta-36.15) ) ) ./ (R .* Ta.^2); % Appendix B
% Pa [Pa] is atmospheric pressure
Pa = 101.3 .*1000 .* ( (293-0.0065.*Elevation) ./293 ).^5.26;
% rou_a [kg/m^3] is the density of air
rou_a = Pa .* Ma ./ (R .* Ta);
% Yita_a [kg/(ms)] is the dynamic viscosity of air
Yita_a = 1.8325 ./ 100000 .* ( 416.16./(Ta+120) ) .* (Ta./296.16).^1.5 ;
% Dv [m^2/s] is the diffusion coefficient for water vapour in air
Dv = 2.11 .* (Ta./273.15).^1.94 .* (Po./Pa) ./100000;
% ka [W/(mK)] is thermal diffusivity of air
ka = 4.1868 ./ 1000 .* (5.96+0.017.*(Ta-273.15));
% Dh [m^2/s] thermal diffusivity of air
Dh = R .* Ta .* ka ./ ( Mw.*lamda.*Gama);
% VPD
%规定时间尺度范围内的饱和水汽压[Pa],在日尺度内的最高与最低气温下：
es_Tmax = 1000 .* 0.6108 .* exp ( (17.27 .* (Tmax - 273.15) ) ./ ( (Tmax - 273.15) + 237.3) );
es_Tmin = 1000 .* 0.6108 .* exp ( (17.27 .* (Tmin - 273.15) ) ./ ( (Tmin - 273.15) + 237.3) );
es = (es_Tmax + es_Tmin) ./ 2;
%由平均相对湿度计算实际水汽压
ea = rhum ./100 .* ((es_Tmax + es_Tmin)./2);
VPD = es-ea;
% 辐射会用到
%！！！原模型的Li是输入变量，中国没有此数据，所以采用FAO的方法，计算Li
Date_str = num2str(Date);
Year_num = str2num(Date_str(:,1:4));
Mon_num = str2num(Date_str(:,5:6));
days = datenum(Year_num,Mon_num,15)-datenum(Year_num,1,1)+1;
%在周期年内地球运转天J日的规格化日-地距离dr和赤道面倾角dt
dr = 1 + 0.033 .* cos(2.*pi.*days./365);
decl = 0.409 .* sin(2*pi.*days/365-1.39);
%太阳高度角由地理纬度f赤道倾角决定
Om = acos(-tan(latd).*tan(decl));
%大气上界太阳辐射Ra，[MJ/(m^2*day)]
Ra = 24*60/pi*0.0820.*dr.*(Om.*sin(latd).*sin(decl) + cos(latd).*cos(decl).*sin(Om));
Ra = Ra.*11.6; %[MJ/(m^2*day)]转化成[W/m^2]
%无云朗空下的地面太阳辐射Rs_o[W/m^2]
Rs_o = (0.75 + 2 * 10^(-5) .* Elevation ).* Ra;
% Sg is global solar irradiance MJ/(m2*day);
sunshinehour_Max = 24.*Om./pi;
Sg = (as+bs.*sunshinehour./sunshinehour_Max).*(Ra./11.6);
Sg = Sg.*11.6; % from （MJ/(m2*day)）to[W/m^2]
% Li incoming long-wave irradiance ！！！原模型这是输入变量，取而代之采用FAO的方法
Li =  Sigma.*Ta.^4 .* (     1 - ( 0.34-0.14.*sqrt(ea./1000) )  .* (1.35*Sg./Rs_o-0.35)     );

%% 3 随蒸发皿变化的参数
% L is the diameter of Class A
D = pan_pars.D; % [m]
L = pan_pars.L; % [m]
% he is the height of rim
he = pan_pars.he; % [m]
% hw is the height of water level
hw = pan_pars.hw; % [m]
% fv需要优化的两个参数
n = pan_pars_n;  % LWH2012中优化出来的参数 ！！！需要优化
q = pan_pars_q; % LWH2012中优化出来的参数 ！！！需要优化
% Beta is the ratio of heat to mass transfer coefficients of the pan
Beta = Dh./Dv .* pan_pars.Beta; % Class A 0.97 外壁面积 0.31 断面面积 1.15 水面面积
% C is the correction factor to account for the shading effect of the bird guard
C = pan_pars.C; %Class A C=1.07; D20 C=1; 601B C=1;
% e_gnd is the emissivity of ground
e_gnd = pan_pars.e_gnd; % LWH2013 这个数可以保持不变
% e_wall is the emissivity of water
e_w = pan_pars.e_w; % LWH2013 这个数可以保持不变
% e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
e_wall =  pan_pars.e_wall; %D20可以使用类似的数据（同一个参考文献），但是601B需要重新确定数值。
% is the refractive index
N = pan_pars.N; % Class A: 1.33
% is the extinction coeeficient
K = pan_pars.K; % Class A: 0
% The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
alpha_0_wall = pan_pars.alpha_0_wall; % Class A : 0.36
% alpha_gnd
alpha_gnd = pan_pars.alpha_gnd;
%% 4 E_pan_A [W/m2] is the aerodynamic component of pan evaporation
% fv
delta_z = L .* ( (rou_a.*L.*n.*U2) ./ (Yita_a) ).^q;
fv = Mw .* Dv ./ ( R .* rou_w .* Ta .* delta_z) ;
% E_pan_A
E_pan_A =  Beta .* Gama .* ( fv .* VPD ) ./ ( Beta .* Gama + s ) ;
E_pan_A(E_pan_A<0) = 0;

%% 5 E_pan_R [W/m2] is the radiative component of pan evaporation
% pan water surface area
Aw = pi .* D.^2 ./ 4;
% Pan water surface area subject to diffuse irradiance
Ad_w = Aw ./ pi .* (  2.*atan(D./he) - he./D.*log(abs( 1+(D./he).^2  ))  );
% Pan wall area subject to diffuse irradiance
Ad_wall = pi.*D.*hw;
% Pan wall area subject to beam irradiance
Ab_wall = D.*hw;
% Outside rim area subject to diffuse irradiance
Ad_rim = pi.*D.*he;
% Outside rim area subject to beam irradiance
Ab_rim = D.*he;
% Inside rim area subject to diffuse irradiance
Ad_rim_pie = Ad_rim.*( D.*log(he./D) + he.*atan(D./he) + 0.5*D.*log(1+D.^2./he.^2) )./(pi.*he);
% Pan bottom area subject to diffuse irradiance
Ad_bot = Aw;
%% 5.1 净长波辐射
% Net long-wave irradiance of the pan water surface [W/m^2]
Ln_w = Li./C .* Ad_w./Aw .* (  1 - (1-e_w) .* (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  )     ) - ...
    e_w .* Sigma .*Ta.^4 .* (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  ) - ...
    e_w .* Sigma .*Ta.^4.*(1-e_wall) .*...
    ( 1 - (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  ) -...
    ( -2.*D.*atan(he./D) - he.*log(abs( 1+(D./he).^2  )) + 4.*D.*atan(0.5.*he./D) + he.*log(abs(1+4.*D.^2./he.^2)) ) ./ (D.*pi)   );
% Net long-wave irradiance of the pan wall
% Ln_wall = e_wall .* Ad_wall ./ Aw .* (  0.5*( (2-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 )  - Sigma.*Ta.^4    ) ;
Ln_rim = e_wall .* Ad_rim ./ Aw .* (  0.5*( (2-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 )  - Sigma.*Ta.^4    ) ;
Ln_rim_pie = 0.5.*Li.*Ad_rim_pie./Aw - e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw .* ( atan(D./he)./pi + 0.5.*D./(he.*pi).*log(abs(1+he.^2./D.^2))   ) - ...
    e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw.*(1-e_wall).* ( -D./(pi.*he).*(-he.*atan(he./D)./D+0.5.*log(abs(1+he.^2./D.^2))) + 2.*D./(pi.*he).*(-0.5.*he.*atan(0.5.*he./D)./D+0.5.*log(abs(1+0.25.*he.^2./D.^2)))  )-...
    e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw.*(1-e_w).* (-he.*atan(D./he) - 0.5.*D.*log(abs(he.^2./D.^2+1)) + 2.*he.*atan(0.5.*D./he) + 0.5.*D.*log(abs(1+4.*he.^2./D.^2))) ./(pi.*he);
% Ln_bot = e_wall .* Ad_bot ./ Aw .* ( (1-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 - Sigma.*Ta.^4 );
Ln = Ln_w  + Ln_rim + Ln_rim_pie ;
%% 5.2 净短波辐射
% fb beam fraction
fb = -0.11 + 1.31.*Sg./Ra;
% Ab_w_pie , alpha_b_w_pie , alpha_d_w , alpha_b_wall_pie , tan_z_pie , alpha_d_wall , alpha_gnd
%Date 是日期信息，201501如2015年1月
oo=1;
for ooo=1:801
    [Ab_w_pie(oo,ooo) , alpha_b_w_pie(oo,ooo) , alpha_d_w(oo,ooo) , alpha_b_wall_pie(oo,ooo) ,...
        tan_z_pie(oo,ooo) , alpha_d_wall(oo,ooo), Ab_rim_pie_pie(oo,ooo)] =...
        S_pars_cal(N,K,D,alpha_0_wall,he,latd(oo,ooo),Date);
    if isnan(alpha_b_w_pie(oo,ooo)) %东北纬度过高导致alpha_b_w_pie算出来NAN
        alpha_b_w_pie(oo,ooo)=0.116;
    end
end
Ab_w_pie=Ab_w_pie(ones(1401,1),:);alpha_b_w_pie=alpha_b_w_pie(ones(1401,1),:);alpha_d_w=alpha_d_w(ones(1401,1),:);
alpha_b_wall_pie=alpha_b_wall_pie(ones(1401,1),:);tan_z_pie=tan_z_pie(ones(1401,1),:);alpha_d_wall=alpha_d_wall(ones(1401,1),:);
Ab_rim_pie_pie=Ab_rim_pie_pie(ones(1401,1),:);

alpha_b_rim_pie = alpha_b_wall_pie;
alpha_d_rim = alpha_d_wall;
alpha_d_bot = alpha_d_wall;
% Net short-wave irradiance of the pan water surface, Sn_w [W/m^2]
if isnan(alpha_b_w_pie) %东北纬度过高导致alpha_b_w_pie算出来NAN
    Sn_w = Sg./(C.*Aw) .* (  fb.*Ab_w_pie.*(1-0.1) + (1-fb).*Ad_w.*(1-alpha_d_w) );
else
    Sn_w = Sg./(C.*Aw) .* (  fb.*Ab_w_pie.*(1-alpha_b_w_pie) + (1-fb).*Ad_w.*(1-alpha_d_w) );
end
% Sn_wall = Sg./Aw .* (     (1-alpha_b_wall_pie).*fb.*tan_z_pie.*Ab_wall    +    Ad_wall.*0.5.*(1-alpha_d_wall).*(alpha_gnd+(1-fb))           );
Sn_rim = Sg./Aw .* (     (1-alpha_b_rim_pie).*fb.*tan_z_pie.*Ab_rim    +    Ad_rim.*0.5.*(1-alpha_d_rim).*(alpha_gnd+(1-fb))           );
Sn_rim_pie = Sg./Aw .* (     fb.*tan_z_pie.*Ab_rim_pie_pie    +    Ad_rim_pie.*0.5.*(1-fb)           );
% Sn_bot = (1-alpha_d_bot).*alpha_gnd.*Sg.*Ad_bot./Aw;
Sn = Sn_w  + Sn_rim + Sn_rim_pie ;
%% 5.3 总净辐射
Rn = Sn + Ln;
E_pan_R = s .* Rn ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
E_pan_R(E_pan_R<0)=0;
Result = E_pan_A + E_pan_R; %m/s
% Rn_lamda_rou是用来画“Figure_E_Rn”文件夹里面的图
Rn_lamda_rou = Rn ./ ( lamda .* rou_w );
% 各个组分分开输出
E_pan_R_w = s .* (Sn_w+Ln_w) ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
E_pan_R_w(E_pan_R_w<0)=0;
% E_pan_R_wall = s .* (Sn_wall+Ln_wall) ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
% E_pan_R_wall(E_pan_R_wall<0)=0;
E_pan_R_rim = s .* (Sn_rim+Ln_rim) ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
E_pan_R_rim(E_pan_R_rim<0)=0;
E_pan_R_rim_pie = s .* (Sn_rim_pie+Ln_rim_pie) ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
E_pan_R_rim_pie(E_pan_R_rim_pie<0)=0;
% E_pan_R_bot = s .* (Sn_bot+Ln_bot) ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
% E_pan_R_bot(E_pan_R_bot<0)=0;
end