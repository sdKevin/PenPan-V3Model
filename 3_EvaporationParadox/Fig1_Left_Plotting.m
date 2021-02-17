function Fig1_Left_Plotting(Epan_Year,YlimRange_Epan,YTickRange_Epan,YTickLabel_Epan,...
    YlimRange_Epan_R,YTickRange_Epan_R,YTickLabel_Epan_R,YlimRange_Epan_A,YTickRange_Epan_A,YTickLabel_Epan_A)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E_pan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Pr
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Epan(1); YlimRange_Epan(1); YlimRange_Epan(2); YlimRange_Epan(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Epan(1); YlimRange_Epan(1); YlimRange_Epan(2); YlimRange_Epan(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan(1) YlimRange_Epan(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Epan
% E_pan m/s to mm/year
Epan_Historical = Epan_Year(1).Epan_Year.E_pan .* 365.*24.*3600.*1000;
Epan_Historical = Epan_Historical - repmat(mean(Epan_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Epan_Historical = nanmean(Epan_Historical)';
c95_Epan_Historical = (std(Epan_Historical)./sqrt(size(Epan_Historical,1))).*1.96; % 95% confidence interval
c95_Epan_Historical = c95_Epan_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Epan_Historical - c95_Epan_Historical; flipud(Ensemble_Mean_Epan_Historical + c95_Epan_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Epan
% E_pan m/s to mm/year
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Epan_Historical = Epan_Year(1).Epan_Year.E_pan .* 365.*24.*3600.*1000;
        Epan_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Epan_Historical = Epan_Year(1).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    end
    Epan_ssp = Epan_Year(i_ssp).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    Epan_ssp = [Epan_Historical(:,end) , Epan_ssp];
    Epan_ssp = Epan_ssp - repmat(mean(Epan_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Epan_ssp(:,i_ssp-1) = nanmean(Epan_ssp)';
    c95_Epan_ssp = (std(Epan_ssp)./sqrt(size(Epan_ssp,1))).*1.96; % 95% confidence interval
    c95_Epan_ssp = c95_Epan_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Epan_ssp(:,i_ssp-1) - c95_Epan_ssp; flipud(Ensemble_Mean_Epan_ssp(:,i_ssp-1) + c95_Epan_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
% Plotting Historical
plot([1850:2014],Ensemble_Mean_Epan_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
% Plotting scenarios
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Epan_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Epan_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Epan_Princeton = Epan_Year(6).Epan_Year.E_pan .* 365.*24.*3600.*1000;
Epan_Princeton = Epan_Princeton - repmat(mean(Epan_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Epan_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Epan_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan(1) YlimRange_Epan(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Epan Anomaly (mm year^-^1)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Epan,'yTick',YTickRange_Epan,'yTickLabel',YTickLabel_Epan,...
    'FontSize',24,'FontName','Arial','TickDir','out','LineWidth',2.5,'XMinorTick','on','YMinorTick','on');
%% Plotting Legend
axes('position',get(gca,'position'),'visible','off')
plot(0,'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot(0,'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot(0,':','Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    end
end
plot(0,':','Color',[238 48 46]./255,'Linewidth',3.5);
legend('historical','ssp126','ssp245','ssp370','ssp585','Princeton',...
    'Location','NorthWest','Color','None','EdgeColor','None','FontSize',24,'FontName','Arial')
set(gca,'visible','off')
clearvars -except Epan_Year YlimRange_Epan YTickRange_Epan YTickLabel_Epan YlimRange_Epan_R YTickRange_Epan_R YTickLabel_Epan_R YlimRange_Epan_A YTickRange_Epan_A YTickLabel_Epan_A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E_pan_R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Pr
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Epan_R(1); YlimRange_Epan_R(1); YlimRange_Epan_R(2); YlimRange_Epan_R(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Epan_R(1); YlimRange_Epan_R(1); YlimRange_Epan_R(2); YlimRange_Epan_R(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan_R(1) YlimRange_Epan_R(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Epan_R
% E_pan_R m/s to mm/year
Epan_R_Historical = Epan_Year(1).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
Epan_R_Historical = Epan_R_Historical - repmat(mean(Epan_R_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Epan_R_Historical = nanmean(Epan_R_Historical)';
c95_Epan_R_Historical = (std(Epan_R_Historical)./sqrt(size(Epan_R_Historical,1))).*1.96; % 95% confidence interval
c95_Epan_R_Historical = c95_Epan_R_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Epan_R_Historical - c95_Epan_R_Historical; flipud(Ensemble_Mean_Epan_R_Historical + c95_Epan_R_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Epan_R
% E_pan_R m/s to mm/year
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Epan_R_Historical = Epan_Year(1).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
        Epan_R_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Epan_R_Historical = Epan_Year(1).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    end
    Epan_R_ssp = Epan_Year(i_ssp).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    Epan_R_ssp = [Epan_R_Historical(:,end) , Epan_R_ssp];
    Epan_R_ssp = Epan_R_ssp - repmat(mean(Epan_R_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Epan_R_ssp(:,i_ssp-1) = nanmean(Epan_R_ssp)';
    c95_Epan_R_ssp = (std(Epan_R_ssp)./sqrt(size(Epan_R_ssp,1))).*1.96; % 95% confidence interval
    c95_Epan_R_ssp = c95_Epan_R_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Epan_R_ssp(:,i_ssp-1) - c95_Epan_R_ssp; flipud(Ensemble_Mean_Epan_R_ssp(:,i_ssp-1) + c95_Epan_R_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Epan_R_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Epan_R_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Epan_R_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Epan_R_Princeton = Epan_Year(6).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
Epan_R_Princeton = Epan_R_Princeton - repmat(mean(Epan_R_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Epan_R_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Epan_R_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan_R(1) YlimRange_Epan_R(2)],'k','LineWidth',1.5);
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Epan_R Anomaly (mm year^-^1)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Epan_R,'yTick',YTickRange_Epan_R,'yTickLabel',YTickLabel_Epan_R,...
    'FontSize',24,'FontName','Arial','TickDir','out','LineWidth',2.5,'XMinorTick','on','YMinorTick','on');
%% Plotting Legend
axes('position',get(gca,'position'),'visible','off')
plot(0,'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot(0,'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot(0,':','Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    end
end
plot(0,':','Color',[238 48 46]./255,'Linewidth',3.5);
legend('historical','ssp126','ssp245','ssp370','ssp585','Princeton',...
    'Location','NorthWest','Color','None','EdgeColor','None','FontSize',24,'FontName','Arial')
set(gca,'visible','off')
clearvars -except Epan_Year YlimRange_Epan YTickRange_Epan YTickLabel_Epan YlimRange_Epan_R YTickRange_Epan_R YTickLabel_Epan_R YlimRange_Epan_A YTickRange_Epan_A YTickLabel_Epan_A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E_pan_A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Pr
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Epan_A(1); YlimRange_Epan_A(1); YlimRange_Epan_A(2); YlimRange_Epan_A(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Epan_A(1); YlimRange_Epan_A(1); YlimRange_Epan_A(2); YlimRange_Epan_A(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan_A(1) YlimRange_Epan_A(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Epan_A
% E_pan_A m/s to mm/year
Epan_A_Historical = Epan_Year(1).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
Epan_A_Historical = Epan_A_Historical - repmat(mean(Epan_A_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Epan_A_Historical = nanmean(Epan_A_Historical)';
c95_Epan_A_Historical = (std(Epan_A_Historical)./sqrt(size(Epan_A_Historical,1))).*1.96; % 95% confidence interval
c95_Epan_A_Historical = c95_Epan_A_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Epan_A_Historical - c95_Epan_A_Historical; flipud(Ensemble_Mean_Epan_A_Historical + c95_Epan_A_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Epan_A
% E_pan_A m/s to mm/year
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Epan_A_Historical = Epan_Year(1).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
        Epan_A_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Epan_A_Historical = Epan_Year(1).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    end
    Epan_A_ssp = Epan_Year(i_ssp).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    Epan_A_ssp = [Epan_A_Historical(:,end) , Epan_A_ssp];
    Epan_A_ssp = Epan_A_ssp - repmat(mean(Epan_A_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Epan_A_ssp(:,i_ssp-1) = nanmean(Epan_A_ssp)';
    c95_Epan_A_ssp = (std(Epan_A_ssp)./sqrt(size(Epan_A_ssp,1))).*1.96; % 95% confidence interval
    c95_Epan_A_ssp = c95_Epan_A_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Epan_A_ssp(:,i_ssp-1) - c95_Epan_A_ssp; flipud(Ensemble_Mean_Epan_A_ssp(:,i_ssp-1) + c95_Epan_A_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Epan_A_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Epan_A_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Epan_A_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Epan_A_Princeton = Epan_Year(6).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
Epan_A_Princeton = Epan_A_Princeton - repmat(mean(Epan_A_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Epan_A_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Epan_A_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Epan_A(1) YlimRange_Epan_A(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Epan_A Anomaly (mm year^-^1)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Epan_A,'yTick',YTickRange_Epan_A,'yTickLabel',YTickLabel_Epan_A,...
    'FontSize',24,'FontName','Arial','TickDir','out','LineWidth',2.5,'XMinorTick','on','YMinorTick','on');
%% Plotting Legend
axes('position',get(gca,'position'),'visible','off')
plot(0,'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot(0,'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot(0,':','Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    end
end
plot(0,':','Color',[238 48 46]./255,'Linewidth',3.5);
legend('historical','ssp126','ssp245','ssp370','ssp585','Princeton',...
    'Location','NorthWest','Color','None','EdgeColor','None','FontSize',24,'FontName','Arial')
set(gca,'visible','off')
clearvars -except Epan_Year YlimRange_Epan YTickRange_Epan YTickLabel_Epan YlimRange_Epan_R YTickRange_Epan_R YTickLabel_Epan_R YlimRange_Epan_A YTickRange_Epan_A YTickLabel_Epan_A

end
