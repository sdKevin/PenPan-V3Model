clc;clear all;
%�в�ƽ������С�Ż��Ľ��
%% 1 DataPrep
[ dem , Tmax_Month , Tmin_Month , U2_Month , SH_Month ,...
    RH_Month , Mask , pars_Efv_601B_n , pars_Efv_601B_q , longd_Mon , latd_Mon , as, bs] = DataPrep();
% 'Latitude(��),','Longitude(��),','Elevation,'��m��
% 'RHU,'(%)��û�о�һ����,'SSD,'(1hr),'MaxT,'(��),'MinT,'(��),'Win'(m/s),as,bs
%% 2 ����������Ĳ��� (601B)
% L is the diameter of Class A
pan_pars.D = 0.618; % [m]
% ע��601B��LӦ������Ȧֱ��
pan_pars.L = pan_pars.D + 0.4; % [m]
% he is the height of rim
pan_pars.he = 0.075; % [m]
% hw is the height of water level
pan_pars.hw = 0.612; % [m]
% Beta is the ratio of heat to mass transfer coefficients of the pan
pan_pars.Beta = 1 + 2*0.618*0.075./(0.25*pi*pan_pars.D^2) + 2*pi*pan_pars.D*0.075./(0.25*pi*pan_pars.D^2); % Class A 0.97 ������ 0.31 ������� 1.15 ˮ�����
% C is the correction factor to account for the shading effect of the bird guard
pan_pars.C = 1; %Class A C=1.07; D20 C=1; 601B C=1;
% e_gnd is the emissivity of ground
pan_pars.e_gnd = 0.90; % LWH2013 ��������Ա��ֲ���
% e_wall is the emissivity of water
pan_pars.e_w = 0.89; % LWH2013 ��������Ա��ֲ���
% e_wall is the emissivity of the pan wallWe assumed wall = 0.82 as quoted for stainless steel coated with zinc oxide (Liebert, 1965, Table I, Substrate Type: B, Coating thickness: 0.05 mm).
pan_pars.e_wall = 0.85; %D20����ʹ�����Ƶ����ݣ�ͬһ���ο����ף���Ϊ���ͬΪZnO��Ĥ������601B��Ҫ����ȷ����ֵ��
% is the refractive index
pan_pars.N = 1.33; % Class A: 1.33
% is the extinction coeeficient
pan_pars.K = 0; % Class A: 0
% The albedo of the pan wall at a solar zenith angle of zero ((Ohman, 1999, Iron galvanised,heavily oxidized))
pan_pars.alpha_0_wall = 0.80; % Class A : 0.36
% alpha_gnd
pan_pars.alpha_gnd = 0.2;
%% 3 ����������
%% 3 ����������
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
            E_pan_R_w(:,:,XuHao),E_pan_R_rim(:,:,XuHao) ,...
            E_pan_R_rim_pie(:,:,XuHao)]...
            = PenPan_V2_Steady_601B(pars_Efv_601B_n(:,:), pars_Efv_601B_q(:,:) , pan_pars, latd_Mon(:,:) , dem(:,:) ,...
            RH_Month(:,:,XuHao) , Tmax_Month(:,:,XuHao) , Tmin_Month(:,:,XuHao) , U2_Month(:,:,XuHao) ,...
            SH_Month(:,:,XuHao) , as(:,:,Month) , bs(:,:,Month) , Date);
        
        Eva_Cal(:,:,XuHao) = Eva_Cal(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        Rn_lamda_rou(:,:,XuHao) = Rn_lamda_rou(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        E_pan_A(:,:,XuHao) = E_pan_A(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        E_pan_R(:,:,XuHao) = E_pan_R(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        E_pan_R_w(:,:,XuHao) = E_pan_R_w(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        E_pan_R_rim(:,:,XuHao) = E_pan_R_rim(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        E_pan_R_rim_pie(:,:,XuHao) = E_pan_R_rim_pie(:,:,XuHao) * Num_Days_Mon(Date) *24 *3600 * 1000;% m/s ת���� mm
        %���������·���Ҫ2��ʱ��
    end
    Year
end
% ע����Ϊrim���ܵķ��������Ȧ/��Ȧ����������������Ҫ��ȥһ�������
Eva_Cal = Eva_Cal - E_pan_R_rim.*0.5 - E_pan_R_rim_pie.*0.5;
E_pan_R_rim = E_pan_R_rim.*0.5;
E_pan_R_rim_pie = E_pan_R_rim_pie.*0.5;
E_pan_R = E_pan_R - E_pan_R_rim.*0.5 - E_pan_R_rim_pie.*0.5;

Eva_Cal(Eva_Cal<0)=0;
Rn_lamda_rou(Rn_lamda_rou<0)=0;
E_pan_A(E_pan_A<0)=0;
E_pan_R(E_pan_R<0)=0;
E_pan_R_w(E_pan_R_w<0)=0;
E_pan_R_rim(E_pan_R_rim<0)=0;
E_pan_R_rim_pie(E_pan_R_rim_pie<0)=0;

save Eva_Cal Eva_Cal
save E_pan_A E_pan_A
save E_pan_R E_pan_R
save E_pan_R_w E_pan_R_w
save E_pan_R_rim E_pan_R_rim
save E_pan_R_rim_pie E_pan_R_rim_pie
%% 4 ��ͼ
Plotting( Eva_Cal )
%% 5 ����Ϊnc����
% ���� Eva_Cal
Value_lon = [70.025:0.05:140.025];
Value_lat = [15.025:0.05:55.025];
Out_Path = ['601B_Epan.nc'];
Object='601B_Epan';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Eva_Cal,Object,days_since,days_end);
% ncdisp(['601B_Epan.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['601B_Epan.nc'], Object);% read variables in nc file
% ���� Epan_wb
Epan_wb = E_pan_A + E_pan_R_w;
Out_Path = ['601B_Epan_wb.nc'];
Object='601B_Epan_wb';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Epan_wb,Object,days_since,days_end);
% ncdisp(['601B_Epan_wb.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['601B_Epan_wb.nc'], Object);% read variables in nc file
% ���� Epan_pw
Epan_pw = E_pan_R - E_pan_R_w;
Out_Path = ['601B_Epan_pw.nc'];
Object='601B_Epan_pw';
days_since = '1960.1';
days_end = '2014.12';
Writ_3Dnc(Out_Path,Value_lon,Value_lat,Epan_pw,Object,days_since,days_end);
% ncdisp(['601B_Epan_pw.nc'],'/', 'full');     % info in nc file
% Variable = ncread(['601B_Epan_pw.nc'], Object);% read variables in nc file