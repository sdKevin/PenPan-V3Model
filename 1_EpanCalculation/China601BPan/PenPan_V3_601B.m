function Epan = PenPan_V3_601B(pan_pars , latd , Elevation , Sg , Ra , Li , U10 , Ta , Sh , Pa)
% Wang, K., Liu, X., Li, Y., Liu, C., & Yang, X. (2018). A generalized evaporation model for Chinese pans.
% Journal of Geophysical Research: Atmospheres, 123,10943иC10966. https://doi.org/10.1029/2018JD028961
Sg(Sg<0)=0; Ra(Ra<0)=0; Li(Li<0)=0; U10(U10<0)=0; Sh(Sh<0)=0; Pa(Pa<0)=0;
% latd бу; Elevation m; Sg W/m2; Li W/m2; U10 m/s; Ta [K]; Pa[Pa]
%% 1 Units and Variables Conversion
latd = pi/180 .* latd; % convert degree to rad
% from U10 to U2 (m/s)
U2 = U10 .* 4.87 ./ (log(67.8.*10-5.42));
clear U10;
%% 2.1 Constants
% Mw [kg/mol] is the molecular mass of water (=0.018 kg/mol).
Mw = 0.018;
% Ma [kg/mol] is the molecular mass of air .
Ma = 0.029;
% R [ J/(mol*K) ] is the ideal gas constant
R = 8.314;
% rou_w [kg/m^3] is the density of liquid water
rou_w = 1000;
% Po [Pa] is the atmospheric pressure at the mean sea level
Po = 101.3 .* 1000;
% Sigma [W/(m^2K^4)] is the Stefan-Boltzmann constant
Sigma = 5.67./100000000;
%% 2.2 Variables for calculating Aerodynamic  component of Epan
% lamda [J/kg] is the latent heat of vaporisation of liquid water
lamda = 2.501 .* 1000000 - 2370 .* ( Ta - 273.15 );
% Gama [Pa/K] psychrometric constant
Gama = 67 - 0.0072 .* Elevation;
% s [Pa/K] is the slope of the saturation vapour pressure versus temperature curve at air temperature (Ta, [K])
s = 611 .* lamda .* Mw .* ( exp( 17.27.*(Ta-273.15)./(Ta-36.15) ) ) ./ (R .* Ta.^2);
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
% VPD [Pa]
ea = Sh .* Pa ./ (0.378 .* Sh + 0.622); % Milly 2016
es = 1000 .* 0.6108 .* exp((17.27.*(Ta-273.15))./Ta);
VPD = es - ea; VPD(VPD<0) = 0;
clear ea es;
%% 3 Pan Parameters
% L is the diameter of a pan
D = pan_pars.D; % [m]
L = pan_pars.L; % [m]
he = pan_pars.he; % [m]
% hw is the height of water level
hw = pan_pars.hw; % [m]
% two parameters in fv
n = (0.7./2) .^ ((0.37-0.0881.*log(U2))./(1-0.0881.*log(0.2))); % n has not adopted the extrapolation of optimized parameters as Shown in Wang et al.2018
q = -0.5; % laminar flow over the pan water surface
% Beta is the ratio of heat to mass transfer coefficients of the pan
Beta = Dh./Dv .* pan_pars.Beta;
% C is the correction factor to account for the shading effect of the bird guard
C = pan_pars.C; % US Class A C=1.07; China D20 C=1; China 601B C=1;
% e_gnd is the emissivity of ground
e_gnd = pan_pars.e_gnd;
% e_wall is the emissivity of water
e_w = pan_pars.e_w;
% e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
e_wall =  pan_pars.e_wall;
% is the refractive index
N = pan_pars.N; % Class A: 1.33
% is the extinction coeeficient
K = pan_pars.K; % Class A: 0
% The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
alpha_0_wall = pan_pars.alpha_0_wall; % Class A : 0.36
% alpha_gnd
alpha_gnd = pan_pars.alpha_gnd;
%% 4 E_pan_A [W/m2] is the aerodynamic component of pan evaporation
% fv
delta_z = L .* ( (rou_a.*L.*n.*U2) ./ (Yita_a) ).^q;
fv = Mw .* Dv ./ ( R .* rou_w .* Ta .* delta_z) ;
% E_pan_A
E_pan_A =  Beta .* Gama .* ( fv .* VPD ) ./ ( Beta .* Gama + s ) ;
E_pan_A(E_pan_A<0) = 0;
clear U2 rou_a Yita_a Dv ka Dh delta_z fv Pa VPD Sh

%% 5 E_pan_R [W/m2] is the radiative component of pan evaporation
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
%% 5.1 Ln [W/m2]: the net long-wave irradiance of the pan.
% Net long-wave irradiance of the pan water surface [W/m^2]
Ln_w = Li./C .* Ad_w./Aw .* (  1 - (1-e_w) .* (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  )     ) - ...
    e_w .* Sigma .*Ta.^4 .* (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  ) - ...
    e_w .* Sigma .*Ta.^4.*(1-e_wall) .*...
    ( 1 - (  2./pi.*atan(D./he) - he.*log(abs( 1+(D./he).^2  ))./(D.*pi)  ) -...
    ( -2.*D.*atan(he./D) - he.*log(abs( 1+(D./he).^2  )) + 4.*D.*atan(0.5.*he./D) + he.*log(abs(1+4.*D.^2./he.^2)) ) ./ (D.*pi)   );
% Net long-wave irradiance of the pan wall, Ln can be smaller than 0
% Ln_wall = e_wall .* Ad_wall ./ Aw .* (  0.5*( (2-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 )  - Sigma.*Ta.^4    ) ;
Ln_rim = e_wall .* Ad_rim ./ Aw .* (  0.5*( (2-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 )  - Sigma.*Ta.^4    ) ;
Ln_rim_pie = 0.5.*Li.*Ad_rim_pie./Aw - e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw .* ( atan(D./he)./pi + 0.5.*D./(he.*pi).*log(abs(1+he.^2./D.^2))   ) - ...
    e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw.*(1-e_wall).* ( -D./(pi.*he).*(-he.*atan(he./D)./D+0.5.*log(abs(1+he.^2./D.^2))) + 2.*D./(pi.*he).*(-0.5.*he.*atan(0.5.*he./D)./D+0.5.*log(abs(1+0.25.*he.^2./D.^2)))  )-...
    e_wall.*Sigma.*Ta.^4.*Ad_rim_pie./Aw.*(1-e_w).* (-he.*atan(D./he) - 0.5.*D.*log(abs(he.^2./D.^2+1)) + 2.*he.*atan(0.5.*D./he) + 0.5.*D.*log(abs(1+4.*he.^2./D.^2))) ./(pi.*he);
% Ln_bot = e_wall .* Ad_bot ./ Aw .* ( (1-e_gnd).*Li + e_gnd.*Sigma.*Ta.^4 - Sigma.*Ta.^4 );
Ln = Ln_w + Ln_rim + Ln_rim_pie;
clear Li Ta
%% 5.2 Sn [W/m2]: the net short-wave irradiance of the pan.
% fb beam fraction
fb = -0.11 + 1.31 .* Sg./Ra;
fb(Ra==0)=0;fb(fb<0)=0;
clear Ra
% Ab_w_pie , alpha_b_w_pie , alpha_d_w , alpha_b_wall_pie , tan_z_pie , alpha_d_wall , alpha_gnd
if exist('Variables_S_pars_cal.mat','file') ==2
    load('Variables_S_pars_cal');
    % Repeat One Year Variables_S_pars_cal to given length
    Ab_w_pie = repmat(Ab_w_pie,[1 1 size(Sg,3)./12]); alpha_b_w_pie = repmat(alpha_b_w_pie,[1 1 size(Sg,3)./12]); alpha_d_w = repmat(alpha_d_w,[1 1 size(Sg,3)./12]);
    alpha_b_wall_pie = repmat(alpha_b_wall_pie,[1 1 size(Sg,3)./12]); tan_z_pie = repmat(tan_z_pie,[1 1 size(Sg,3)./12]); alpha_d_wall = repmat(alpha_d_wall,[1 1 size(Sg,3)./12]);
    Ab_rim_pie_pie = repmat(Ab_rim_pie_pie,[1 1 size(Sg,3)./12]);
else
    for i_mon = 1 : 12
        % Date represents Jan-1850 as 185001
        Year = ceil(i_mon./12) + 1849;
        Month = mod(i_mon,12);
        if Month==0
            Month = 12;
        end
        if Month<10
            Date = [num2str(Year) '0' num2str(Month)];
            Date = str2num(Date);
        else
            Date = [num2str(Year) num2str(Month)];
            Date = str2num(Date);
        end
        % Calculate Variables_S_pars_cal
        oo = 1;
        for ooo = 1 : size(latd,2)
            [Ab_w_pie(oo,ooo,i_mon) , alpha_b_w_pie(oo,ooo,i_mon) , alpha_d_w(oo,ooo,i_mon) , alpha_b_wall_pie(oo,ooo,i_mon) ,...
                tan_z_pie(oo,ooo,i_mon) , alpha_d_wall(oo,ooo,i_mon), Ab_rim_pie_pie(oo,ooo,i_mon)] =...
                S_pars_cal(N , K , D , alpha_0_wall , he , latd(oo,ooo) , Date);
        end
        % For a latitude, all these variables are the same
        if i_mon==1
            Ab_w_pie=Ab_w_pie(ones(720,1),:,i_mon);alpha_b_w_pie=alpha_b_w_pie(ones(720,1),:,i_mon);alpha_d_w=alpha_d_w(ones(720,1),:,i_mon);
            alpha_b_wall_pie=alpha_b_wall_pie(ones(720,1),:,i_mon);tan_z_pie=tan_z_pie(ones(720,1),:,i_mon);alpha_d_wall=alpha_d_wall(ones(720,1),:,i_mon);
            Ab_rim_pie_pie=Ab_rim_pie_pie(ones(720,1),:,i_mon);
        else
            Ab_w_pie(:,:,i_mon)=Ab_w_pie(ones(720,1),:,i_mon);alpha_b_w_pie(:,:,i_mon)=alpha_b_w_pie(ones(720,1),:,i_mon);alpha_d_w(:,:,i_mon)=alpha_d_w(ones(720,1),:,i_mon);
            alpha_b_wall_pie(:,:,i_mon)=alpha_b_wall_pie(ones(720,1),:,i_mon);tan_z_pie(:,:,i_mon)=tan_z_pie(ones(720,1),:,i_mon);alpha_d_wall(:,:,i_mon)=alpha_d_wall(ones(720,1),:,i_mon);
            Ab_rim_pie_pie(:,:,i_mon)=Ab_rim_pie_pie(ones(720,1),:,i_mon);
        end
    end
    
    save('Variables_S_pars_cal' , 'Ab_w_pie' , 'alpha_b_w_pie' , 'alpha_d_w',...
        'alpha_b_wall_pie' , 'tan_z_pie' , 'alpha_d_wall' , 'Ab_rim_pie_pie');
    % Repeat One Year Variables_S_pars_cal to given length
    Ab_w_pie = repmat(Ab_w_pie,[1 1 size(Sg,3)./12]); alpha_b_w_pie = repmat(alpha_b_w_pie,[1 1 size(Sg,3)./12]); alpha_d_w = repmat(alpha_d_w,[1 1 size(Sg,3)./12]);
    alpha_b_wall_pie = repmat(alpha_b_wall_pie,[1 1 size(Sg,3)./12]); tan_z_pie = repmat(tan_z_pie,[1 1 size(Sg,3)./12]); alpha_d_wall = repmat(alpha_d_wall,[1 1 size(Sg,3)./12]);
    Ab_rim_pie_pie = repmat(Ab_rim_pie_pie,[1 1 size(Sg,3)./12]);
end
alpha_b_rim_pie = alpha_b_wall_pie;
alpha_d_rim = alpha_d_wall;
alpha_d_bot = alpha_d_wall;
% Net short-wave irradiance of the pan water surface, Sn_w [W/m^2]
Sn_w = Sg./(C.*Aw) .* (  fb.*Ab_w_pie.*(1-alpha_b_w_pie) + (1-fb).*Ad_w.*(1-alpha_d_w) );
Sn_w(Sn_w<0)=0;
% Sn_wall = Sg./Aw .* (     (1-alpha_b_wall_pie).*fb.*tan_z_pie.*Ab_wall    +    Ad_wall.*0.5.*(1-alpha_d_wall).*(alpha_gnd+(1-fb))           );
% Sn_wall(Sn_wall<0)=0;
Sn_rim = Sg./Aw .* (     (1-alpha_b_rim_pie).*fb.*tan_z_pie.*Ab_rim    +    Ad_rim.*0.5.*(1-alpha_d_rim).*(alpha_gnd+(1-fb))           );
Sn_rim(Sn_rim<0)=0;
Sn_rim_pie = Sg./Aw .* (     fb.*tan_z_pie.*Ab_rim_pie_pie    +    Ad_rim_pie.*0.5.*(1-fb)           );
Sn_rim_pie(Sn_rim_pie<0)=0;
% Sn_bot = (1-alpha_d_bot).*alpha_gnd.*Sg.*Ad_bot./Aw;
% Sn_bot(Sn_bot<0)=0;
Sn = Sn_w + Sn_rim + Sn_rim_pie;
clear fb Sg Ab_w_pie alpha_b_w_pie alpha_d_w alpha_b_wall_pie tan_z_pie alpha_d_wall Ab_rim_pie_pie alpha_b_rim_pie alpha_d_rim alpha_d_bot
%% 5.3 Rn [W/m2]: the net irradiance of the pan.
Rn = Sn + Ln;
E_pan_R = s .* Rn ./ ( ( s + Beta .* Gama ) .* lamda .* rou_w );
E_pan_R(E_pan_R<0)=0;
E_pan = E_pan_A + E_pan_R; %m/s
% Different Components
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

%% Output
Epan.E_pan = E_pan - E_pan_R_rim .* 0.5 - E_pan_R_rim_pie .* 0.5; % Due to water in double-ring and measuring the inner one.
Epan.E_pan_R = E_pan_R - E_pan_R_rim .* 0.5 - E_pan_R_rim_pie .* 0.5;
Epan.E_pan_A = E_pan_A;
Epan.E_pan_R_w = E_pan_R_w;
% Epan.E_pan_R_wall = E_pan_R_wall;
Epan.E_pan_R_rim = E_pan_R_rim .* 0.5;
Epan.E_pan_R_rim_pie = E_pan_R_rim_pie .* 0.5;
% Epan.E_pan_R_bot = E_pan_R_bot;
end