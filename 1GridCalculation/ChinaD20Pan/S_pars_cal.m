function [Ab_w_pie , alpha_b_w_pie , alpha_d_w , alpha_b_wall_pie , tan_z_pie , alpha_d_wall , Ab_rim_pie_pie] =...
    S_pars_cal(N,K,D,alpha_0_wall,he,latd,Date)
% decl 某一天的declination
% !!! 用integral进行积分不会报错，因为可以被解出来；不需要用quad进行数值积分，而且quad还会在间断点处取值出现问题；
%% 1 Ab_w_pie
% 计算该月份有几天
Num_days = Num_Days_Mon(Date);
for iii = 1:Num_days
    %根据月份和第几天算出decl
    Date_str = num2str(Date);
    Year_num = str2num(Date_str(:,1:4));
    Mon_num = str2num(Date_str(:,5:6));
    days = datenum(Year_num,Mon_num,iii)-datenum(Year_num,1,1)+1;
    %Roderick 1999 % GGG = (days-1)*2*pi/365;
    %Roderick 1999 % decl = 0.006918-0.399912*cos(GGG)+0.0770257*sin(GGG)-0.006758*cos(2*GGG)+0.000907*sin(2*GGG)-0.002697*cos(3*GGG)+0.00148*sin(3*GGG);
    decl =  0.409 .* sin(2*pi.*days/365-1.39);
    %半天长计算
    day_len = acos(-tan(latd).*tan(decl)).*2./(15.*pi./180); %the number of daylight hours (day_len,hour) for a particular day using the hour angle at sunset  
    day_rad = day_len.*2.*pi./24; %from hour to rad
    %积分下限
    omega_s_ = -0.5.*day_rad;
    %积分上限
    omega_s = 0.5.*day_rad;
    % 获取被积函数的显式表达形式！
    %         syms hourang_rad decl latd D he% 积分变量
    %         cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(hourang_rad);
    %         z = acos(cosz);
    %         tanz = tan(z);
    %         Ab_w = 0.5.*( D.^2.*acos(he.*tanz./D) - he.*tanz.*sqrt(D.^2-he.^2.*tanz.^2));
    %         Inte1 =  Ab_w.*cosz.*86400./(2.*pi);
    Inte1_inline = @(hourang_rad) ( ((1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1/2)./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) < (D./he)).*...
        ((43200.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).*((D.^2.*acos((he.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(D.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./2 - (he.*(D.^2 + (he.^2.*((sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2 - 1))./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2).*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(2.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./pi);
    Inte1_result(iii) = integral(Inte1_inline,omega_s_,omega_s)./day_rad;
    
    % 获取被积函数的显式表达形式！
    %     cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(hourang_rad);
    %     Inte2 = cosz.*86400./(2.*pi);
    Inte2_inline =  @(hourang_rad) ((86400.*sin(decl).*sin(latd) + 86400.*cos(decl).*cos(hourang_rad).*cos(latd))./(2.*pi));
    Inte2_result(iii) = integral(Inte2_inline,omega_s_,omega_s)./day_rad;
end
Ab_w_pie = sum(Inte1_result)./sum(Inte2_result);
clear iii Inte1 Inte1_inline  Inte1_result Inte2 Inte2_inline Inte2_result

%% 2 alpha_b_w_pie
alpha_0 = ( (N-1).^2 + K.^2 ) ./ ( (N+1).^2 + K.^2 )  ;
for iii = 1:Num_days
    % 根据月份和第几天算出decl
    Date_str = num2str(Date);
    Year_num = str2num(Date_str(:,1:4));
    Mon_num = str2num(Date_str(:,5:6));
    days = datenum(Year_num,Mon_num,iii)-datenum(Year_num,1,1)+1;
    decl =  0.409 .* sin(2*pi.*days/365-1.39);
    % 半天长计算
    day_len = acos(-tan(latd).*tan(decl)).*2./(15.*pi./180);
    day_rad = day_len.*2.*pi./24;
    % 积分下限
    omega_s_ = -0.5.*day_rad;
    % 积分上限
    omega_s = 0.5.*day_rad;
    % 获取被积函数的显式表达形式！
    %         syms hourang_rad decl latd alpha_0 D he% 积分变量
    %         cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(hourang_rad);
    %         z = acos(cosz);
    %         Ab_w = 0.5.*( D.^2.*acos(he.*tan(z)./D) - he.*tan(z).*sqrt(D.^2-he.^2.*tan(z).^2));
    %         alpha_b_w = alpha_0+(1-alpha_0).*(1-cosz).^5;
    %         Inte1 = alpha_b_w.*Ab_w.*cosz.*86400./(2.*pi);
    Inte1_inline = @(hourang_rad) ( ((1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1/2)./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) < (D./he)).*...
        ((43200.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).*(alpha_0 + (alpha_0 - 1).*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd) - 1).^5).*((D.^2.*acos((he.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(D.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./2 - (he.*(D.^2 + (he.^2.*((sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2 - 1))./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2).*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(2.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./pi);
    Inte1_result(iii) = integral(Inte1_inline,omega_s_,omega_s)./day_rad;
    Inte2_inline = @(hourang_rad) ( ((1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1/2)./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) < (D./he)).*...
        ((43200.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).*((D.^2.*acos((he.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(D.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./2 - (he.*(D.^2 + (he.^2.*((sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2 - 1))./(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2).*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./(2.*(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))))./pi);
    Inte2_result(iii) = integral(Inte2_inline,omega_s_,omega_s)./day_rad;
end
alpha_b_w_pie = sum(Inte1_result)./sum(Inte2_result);
clear iii Inte1 Inte1_inline Inte1_result Inte2 Inte2_inline Inte2_result

%% 3 alpha_d_w
% LWH原文公式(C.2)有问题，需要说明进行修改
alpha_0 = ( (N-1).^2 + K.^2 ) ./ ( (N+1).^2 + K.^2 )  ; %equation (23) in LWH 2013
% syms x D he alpha_0
% theta_1 = atan((D-x)./he);
% theta_2 = atan(x./he);
% AAA = alpha_0 + (1-alpha_0) ./ (160.*theta_1) .* ( 1260.*theta_1 - 2100.*sin(theta_1) + ...
%     600.*sin(2.*theta_1) - 150.*sin(3.*theta_1) + 25.*sin(4.*theta_1) - 2.*sin(5.*theta_1) );
% BBB = alpha_0 + (1-alpha_0) ./ (160.*theta_2) .* ( 1260.*theta_2 - 2100.*sin(theta_2) + ...
%     600.*sin(2.*theta_2) - 150.*sin(3.*theta_2) + 25.*sin(4.*theta_2) - 2.*sin(5.*theta_2) );
% Integrand = (theta_1.*AAA+theta_2.*BBB)./(theta_1+theta_2);
Integrand_inline = @(x) (atan(x./he).*(alpha_0 - ((alpha_0 - 1).*(1260.*atan(x./he) + 600.*sin(2.*atan(x./he)) - 150.*sin(3.*atan(x./he)) + 25.*sin(4.*atan(x./he)) - 2.*sin(5.*atan(x./he)) - (2100.*x)./(he.*(x.^2./he.^2 + 1).^(1./2))))./(160.*atan(x./he))) + atan((D - x)./he).*(alpha_0 - ((alpha_0 - 1).*(1260.*atan((D - x)./he) + 600.*sin(2.*atan((D - x)./he)) - 150.*sin(3.*atan((D - x)./he)) + 25.*sin(4.*atan((D - x)./he)) - 2.*sin(5.*atan((D - x)./he)) - (2100.*(D - x))./(he.*((D - x).^2./he.^2 + 1).^(1./2))))./(160.*atan((D - x)./he))))./(atan(x./he) + atan((D - x)./he));
alpha_d_w = integral(Integrand_inline,0,D)./D;

%% 4 alpha_b_wall_pie
for iii = 1:Num_days
    %根据月份和第几天算出decl
    Date_str = num2str(Date);
    Year_num = str2num(Date_str(:,1:4));
    Mon_num = str2num(Date_str(:,5:6));
    days = datenum(Year_num,Mon_num,iii)-datenum(Year_num,1,1)+1;
    decl =  0.409 .* sin(2*pi.*days/365-1.39);
    %半天长计算
    day_len = acos(-tan(latd).*tan(decl)).*2./(15.*pi./180);
    day_rad = day_len.*2.*pi./24;
    %积分下限
    omega_s_ = -0.5.*day_rad;
    %积分上限
    omega_s = 0.5.*day_rad;
    %被积函数
    %     syms decl latd hourang_rad alpha_0_wall % 积分变量
    %     cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(hourang_rad);
    %     z = acos(cosz);
    %     sinz = sin(z);
    %     alpha_b_wall = alpha_0_wall + (1-alpha_0_wall).*(1260.*z+2100.*cos(z)-600.*sin(2.*z)-150.*cos(3.*z)+25.*sin(4.*z)+2.*cos(5.*z)-1952)./(160.*z);
    %     Inte1 = alpha_b_wall.*sinz.*86400./(2.*pi);
    Inte1_inline = @(hourang_rad) (43200.*(alpha_0_wall - ((alpha_0_wall - 1).*(2.*cos(5.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) - 150.*cos(3.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) - 600.*sin(2.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) + 25.*sin(4.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd))) + 1260.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)) + 2100.*sin(decl).*sin(latd) + 2100.*cos(decl).*cos(hourang_rad).*cos(latd) - 1952))./(160.*acos(sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)))).*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./pi;
    Inte1_result(iii) = integral(Inte1_inline,omega_s_,omega_s)./day_rad;
    %     Inte2 = sinz.*86400./(2.*pi);
    Inte2_inline = @(hourang_rad) (43200.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./pi;
    Inte2_result(iii) = integral(Inte2_inline,omega_s_,omega_s)./day_rad;
end
alpha_b_wall_pie = sum(Inte1_result)./sum(Inte2_result);

%% 5 tan_z_pie
for iii = 1:Num_days
    %根据月份和第几天算出decl
    Date_str = num2str(Date);
    Year_num = str2num(Date_str(:,1:4));
    Mon_num = str2num(Date_str(:,5:6));
    days = datenum(Year_num,Mon_num,iii)-datenum(Year_num,1,1)+1;
    decl =  0.409 .* sin(2*pi.*days/365-1.39);
    %半天长计算
    day_len = acos(-tan(latd).*tan(decl)).*2./(15.*pi./180);
    day_rad = day_len.*2.*pi./24;
    %积分下限
    omega_s_ = -0.5.*day_rad;
    %积分上限
    omega_s = 0.5.*day_rad;
    %被积函数
    %     syms decl latd hourang_rad % 积分变量
    %     cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(hourang_rad);
    %     z = acos(cosz);
    %     sinz = sin(z);
    %     Inte1 = sinz.*86400./(2.*pi);
    Inte1_inline = @(hourang_rad) (43200.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./pi;
    Inte1_result(iii) = integral(Inte1_inline,omega_s_,omega_s)./day_rad;
    %     Inte2 = cosz.*86400./(2.*pi);
    Inte2_inline = @(hourang_rad) (86400.*sin(decl).*sin(latd) + 86400.*cos(decl).*cos(hourang_rad).*cos(latd))./(2.*pi);
    Inte2_result(iii) = integral(Inte2_inline,omega_s_,omega_s)./day_rad;
end
tan_z_pie = sum(Inte1_result)./sum(Inte2_result);

%% 6 alpha_d_wall
alpha_d_wall = alpha_0_wall + (1-alpha_0_wall).*(315*pi-976)./(40*pi);

%% Ab_rim_pie_pie
% 计算该月份有几天
Num_days = Num_Days_Mon(Date);
for iii = 1:Num_days
    %根据月份和第几天算出decl
    Date_str = num2str(Date);
    Year_num = str2num(Date_str(:,1:4));
    Mon_num = str2num(Date_str(:,5:6));
    days = datenum(Year_num,Mon_num,iii)-datenum(Year_num,1,1)+1;
    %Roderick 1999 % GGG = (days-1)*2*pi/365;
    %Roderick 1999 % decl = 0.006918-0.399912*cos(GGG)+0.0770257*sin(GGG)-0.006758*cos(2*GGG)+0.000907*sin(2*GGG)-0.002697*cos(3*GGG)+0.00148*sin(3*GGG);
    decl =  0.409 .* sin(2*pi.*days/365-1.39);
    %半天长计算
    day_len = acos(-tan(latd).*tan(decl)).*2./(15.*pi./180);
    day_rad = day_len.*2.*pi./24;
    %积分下限
    omega_s_ = -0.5.*day_rad;
    %积分上限
    omega_s = 0.5.*day_rad;
    for iiii = 1:1000
        iiii_w = omega_s_ + day_rad*iiii./1000;
        cosz = sin(decl).*sin(latd) + cos(decl).*cos(latd).*cos(iiii_w);
        z = acos(cosz);
        while z>(0.5*pi)
            z = 0.5*pi;
        end
        while z<0
            z = 0;
        end
        sinz = sin(z);
        tanz = tan(z);
        if tanz<(D/he)
            Ab_rim_pie = 0.25.*pi.*D.^2./tanz + 0.5*he.*sqrt(D.^2-tanz.^2.*he.^2) - 0.5.*D.^2.*asin(sqrt(1-tanz.^2.*he.^2./D.^2))./tanz;
        else
            Ab_rim_pie = 0.25.*pi.*D.^2./tanz;
        end
        Inte1(iiii) =  Ab_rim_pie.*sinz.*86400./(2.*pi);
        Inte2(iiii) =  sinz.*86400./(2.*pi);
    end
    Inte1_result(iii) = sum(Inte1)./day_rad;
    Inte2_result(iii) = sum(Inte2)./day_rad;
    %     Inte2_inline = @(hourang_rad) (43200.*(1 - (sin(decl).*sin(latd) + cos(decl).*cos(hourang_rad).*cos(latd)).^2).^(1./2))./pi;
    %     Inte2_result(iii) = integral(Inte2_inline,omega_s_,omega_s)./day_rad;
end
Ab_rim_pie_pie = sum(Inte1_result)./sum(Inte2_result);
clear iii Inte1 Inte1_inline  Inte1_result Inte2 Inte2_inline Inte2_result
end