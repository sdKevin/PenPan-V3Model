function Fig2_Left_Plotting(Met_Year)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Sg
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Sg = [-6,6]; YTickRange_Sg = [-6:2:6]; YTickLabel_Sg = {'';'-4';'';'0';'';'4';''};
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Sg(1); YlimRange_Sg(1); YlimRange_Sg(2); YlimRange_Sg(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Sg(1); YlimRange_Sg(1); YlimRange_Sg(2); YlimRange_Sg(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Sg(1) YlimRange_Sg(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Sg
% Sg W/m2
Sg_Historical = Met_Year(1).Met_Year.Sg;
Sg_Historical = Sg_Historical - repmat(mean(Sg_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Sg_Historical = nanmean(Sg_Historical)';
c95_Sg_Historical = (std(Sg_Historical)./sqrt(size(Sg_Historical,1))).*1.96; % 95% confidence interval
c95_Sg_Historical = c95_Sg_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Sg_Historical - c95_Sg_Historical; flipud(Ensemble_Mean_Sg_Historical + c95_Sg_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Sg
% Sg W/m2
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Sg_Historical = Met_Year(1).Met_Year.Sg;
        Sg_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Sg_Historical = Met_Year(1).Met_Year.Sg;
    end
    Sg_ssp = Met_Year(i_ssp).Met_Year.Sg;
    Sg_ssp = [Sg_Historical(:,end) , Sg_ssp];
    Sg_ssp = Sg_ssp - repmat(mean(Sg_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Sg_ssp(:,i_ssp-1) = nanmean(Sg_ssp)';
    c95_Sg_ssp = (std(Sg_ssp)./sqrt(size(Sg_ssp,1))).*1.96; % 95% confidence interval
    c95_Sg_ssp = c95_Sg_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Sg_ssp(:,i_ssp-1) - c95_Sg_ssp; flipud(Ensemble_Mean_Sg_ssp(:,i_ssp-1) + c95_Sg_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Sg_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Sg_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Sg_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Sg_Princeton = Met_Year(6).Met_Year.Sg;
Sg_Princeton = Sg_Princeton - repmat(mean(Sg_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Sg_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Sg_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Sg(1) YlimRange_Sg(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Sg Anomaly (W m^-^2)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Sg,'yTick',YTickRange_Sg,'yTickLabel',YTickLabel_Sg,...
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
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Li %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Li
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Li = [-10,50]; YTickRange_Li = [-10:10:50]; YTickLabel_Li = {'';'0';'';'20';'';'40';''};
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Li(1); YlimRange_Li(1); YlimRange_Li(2); YlimRange_Li(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Li(1); YlimRange_Li(1); YlimRange_Li(2); YlimRange_Li(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Li(1) YlimRange_Li(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Li
% Li W/m2
Li_Historical = Met_Year(1).Met_Year.Li;
Li_Historical = Li_Historical - repmat(mean(Li_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Li_Historical = nanmean(Li_Historical)';
c95_Li_Historical = (std(Li_Historical)./sqrt(size(Li_Historical,1))).*1.96; % 95% confidence interval
c95_Li_Historical = c95_Li_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Li_Historical - c95_Li_Historical; flipud(Ensemble_Mean_Li_Historical + c95_Li_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Li
% Li W/m2
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Li_Historical = Met_Year(1).Met_Year.Li;
        Li_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Li_Historical = Met_Year(1).Met_Year.Li;
    end
    Li_ssp = Met_Year(i_ssp).Met_Year.Li;
    Li_ssp = [Li_Historical(:,end) , Li_ssp];
    Li_ssp = Li_ssp - repmat(mean(Li_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Li_ssp(:,i_ssp-1) = nanmean(Li_ssp)';
    c95_Li_ssp = (std(Li_ssp)./sqrt(size(Li_ssp,1))).*1.96; % 95% confidence interval
    c95_Li_ssp = c95_Li_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Li_ssp(:,i_ssp-1) - c95_Li_ssp; flipud(Ensemble_Mean_Li_ssp(:,i_ssp-1) + c95_Li_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Li_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Li_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Li_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Li_Princeton = Met_Year(6).Met_Year.Li;
Li_Princeton = Li_Princeton - repmat(mean(Li_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Li_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Li_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Li(1) YlimRange_Li(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Li Anomaly (W m^-^2)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Li,'yTick',YTickRange_Li,'yTickLabel',YTickLabel_Li,...
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
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% U10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% U10
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_U10 = [-0.08,0.06]; YTickRange_U10 = [-0.08:0.04:0.06]; YTickLabel_U10 = {'-0.08';'-0.04';'0';'0.04'};
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_U10(1); YlimRange_U10(1); YlimRange_U10(2); YlimRange_U10(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_U10(1); YlimRange_U10(1); YlimRange_U10(2); YlimRange_U10(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_U10(1) YlimRange_U10(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical U10
% U10 W/m2
U10_Historical = Met_Year(1).Met_Year.U10;
U10_Historical = U10_Historical - repmat(mean(U10_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_U10_Historical = nanmean(U10_Historical)';
c95_U10_Historical = (std(U10_Historical)./sqrt(size(U10_Historical,1))).*1.96; % 95% confidence interval
c95_U10_Historical = c95_U10_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_U10_Historical - c95_U10_Historical; flipud(Ensemble_Mean_U10_Historical + c95_U10_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario U10
% U10 W/m2
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        U10_Historical = Met_Year(1).Met_Year.U10;
        U10_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        U10_Historical = Met_Year(1).Met_Year.U10;
    end
    U10_ssp = Met_Year(i_ssp).Met_Year.U10;
    U10_ssp = [U10_Historical(:,end) , U10_ssp];
    U10_ssp = U10_ssp - repmat(mean(U10_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_U10_ssp(:,i_ssp-1) = nanmean(U10_ssp)';
    c95_U10_ssp = (std(U10_ssp)./sqrt(size(U10_ssp,1))).*1.96; % 95% confidence interval
    c95_U10_ssp = c95_U10_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_U10_ssp(:,i_ssp-1) - c95_U10_ssp; flipud(Ensemble_Mean_U10_ssp(:,i_ssp-1) + c95_U10_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_U10_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_U10_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_U10_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
U10_Princeton = Met_Year(6).Met_Year.U10;
U10_Princeton = U10_Princeton - repmat(mean(U10_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],U10_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],U10_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_U10(1) YlimRange_U10(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('U10 Anomaly (m s^-^1)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_U10,'yTick',YTickRange_U10,'yTickLabel',YTickLabel_U10,...
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
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Ta
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Ta = [-2,8]; YTickRange_Ta = [-2:2:8]; YTickLabel_Ta = {'';'0';'';'4';'';'8'};
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Ta(1); YlimRange_Ta(1); YlimRange_Ta(2); YlimRange_Ta(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Ta(1); YlimRange_Ta(1); YlimRange_Ta(2); YlimRange_Ta(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Ta(1) YlimRange_Ta(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Ta
% Ta W/m2
Ta_Historical = Met_Year(1).Met_Year.Ta;
Ta_Historical = Ta_Historical - repmat(mean(Ta_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Ta_Historical = nanmean(Ta_Historical)';
c95_Ta_Historical = (std(Ta_Historical)./sqrt(size(Ta_Historical,1))).*1.96; % 95% confidence interval
c95_Ta_Historical = c95_Ta_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Ta_Historical - c95_Ta_Historical; flipud(Ensemble_Mean_Ta_Historical + c95_Ta_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Ta
% Ta W/m2
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Ta_Historical = Met_Year(1).Met_Year.Ta;
        Ta_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Ta_Historical = Met_Year(1).Met_Year.Ta;
    end
    Ta_ssp = Met_Year(i_ssp).Met_Year.Ta;
    Ta_ssp = [Ta_Historical(:,end) , Ta_ssp];
    Ta_ssp = Ta_ssp - repmat(mean(Ta_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Ta_ssp(:,i_ssp-1) = nanmean(Ta_ssp)';
    c95_Ta_ssp = (std(Ta_ssp)./sqrt(size(Ta_ssp,1))).*1.96; % 95% confidence interval
    c95_Ta_ssp = c95_Ta_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Ta_ssp(:,i_ssp-1) - c95_Ta_ssp; flipud(Ensemble_Mean_Ta_ssp(:,i_ssp-1) + c95_Ta_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Ta_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Ta_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Ta_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Ta_Princeton = Met_Year(6).Met_Year.Ta;
Ta_Princeton = Ta_Princeton - repmat(mean(Ta_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Ta_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Ta_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Ta(1) YlimRange_Ta(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Ta Anomaly (K)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Ta,'yTick',YTickRange_Ta,'yTickLabel',YTickLabel_Ta,...
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
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%% Setting Color and axis property
RGB_Historical_Shade = [205,205,205]; RGB_Historical_Line = [23,23,23];
% Sh
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Sh = [-0.5,3.5]./1000; YTickRange_Sh = [-0.5:0.5:3.5]./1000; YTickLabel_Sh = {'';'0';'';'1';'';'2';'';'3';''};
%% Three time windows
%  Contemporary: 1948-2014
fill([1948;2014;2014;1948],...
    [YlimRange_Sh(1); YlimRange_Sh(1); YlimRange_Sh(2); YlimRange_Sh(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
% Far Future: 2070-2099
fill([2070;2099;2099;2070],...
    [YlimRange_Sh(1); YlimRange_Sh(1); YlimRange_Sh(2); YlimRange_Sh(2)],...
    [240,240,242]./255,'EdgeAlpha',0,'FaceAlpha',0.9);
%% Plotting Shade Area
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Sh(1) YlimRange_Sh(2)],'k','LineWidth',1.5); hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% Historical Sh
% Sh
Sh_Historical = Met_Year(1).Met_Year.Sh;
Sh_Historical = Sh_Historical - repmat(mean(Sh_Historical(:,99:165),2),1,165); % Change to Anomaly
Ensemble_Mean_Sh_Historical = nanmean(Sh_Historical)';
c95_Sh_Historical = (std(Sh_Historical)./sqrt(size(Sh_Historical,1))).*1.96; % 95% confidence interval
c95_Sh_Historical = c95_Sh_Historical';
h1 = fill([[1850:2014]';flipud([1850:2014]')],...
    [Ensemble_Mean_Sh_Historical - c95_Sh_Historical; flipud(Ensemble_Mean_Sh_Historical + c95_Sh_Historical)],...
    RGB_Historical_Shade./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;

% Scenario Sh
% Sh
for i_ssp = [5,4,3,2]
    if i_ssp == 4
        Sh_Historical = Met_Year(1).Met_Year.Sh;
        Sh_Historical(16,:) = []; % HadGEM3-GC31-LL
    else
        Sh_Historical = Met_Year(1).Met_Year.Sh;
    end
    Sh_ssp = Met_Year(i_ssp).Met_Year.Sh;
    Sh_ssp = [Sh_Historical(:,end) , Sh_ssp];
    Sh_ssp = Sh_ssp - repmat(mean(Sh_Historical(:,99:165),2), 1 ,87); % Change to Anomaly
    Ensemble_Mean_Sh_ssp(:,i_ssp-1) = nanmean(Sh_ssp)';
    c95_Sh_ssp = (std(Sh_ssp)./sqrt(size(Sh_ssp,1))).*1.96; % 95% confidence interval
    c95_Sh_ssp = c95_Sh_ssp';
    if i_ssp == 5|| i_ssp == 2
        h_ssp = fill([[2014:2100]';flipud([2014:2100]')],...
            [Ensemble_Mean_Sh_ssp(:,i_ssp-1) - c95_Sh_ssp; flipud(Ensemble_Mean_Sh_ssp(:,i_ssp-1) + c95_Sh_ssp)],...
            RGB_ssp_Shade(i_ssp-1,:)./255,'EdgeAlpha',0,'FaceAlpha',0.9); hold on;
    end
end
%% Plotting Ensemble Mean
plot([1850:2014],Ensemble_Mean_Sh_Historical,'-',...
    'Color',RGB_Historical_Line./255,'Linewidth',3.5);hold on;
for i_ssp = [1,2,3,4]
    if i_ssp == 1|| i_ssp ==4
        plot([2014:2100],Ensemble_Mean_Sh_ssp(:,i_ssp),'-',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',3.5);
    else
        plot([2014:2100],Ensemble_Mean_Sh_ssp(:,i_ssp),':',...
            'Color',RGB_ssp_Line(i_ssp,:)./255,'Linewidth',2.5);
    end
end
% Plotting Princeton
Sh_Princeton = Met_Year(6).Met_Year.Sh;
Sh_Princeton = Sh_Princeton - repmat(mean(Sh_Princeton(:)),1,67); % Change to Anomaly
plot([1948:2014],Sh_Princeton,'Color',[238 48 46]./255,'Linewidth',3.5);
plot([1948:2014],Sh_Princeton,'--','Color',[1 1 1],'Linewidth',3.5);
% Plot y=0 and x=2014
plot([2014 2014],[YlimRange_Sh(1) YlimRange_Sh(2)],'k','LineWidth',1.5);
hold on;
plot([1948 2100],[0 0],'Color',[189,188,188]./255,'LineWidth',3)
% setting axis
ylabel('Sh Anomaly (¡Á10^-^3)');
set(gca,'xlim',[1948,2100],'ylim',YlimRange_Sh,'yTick',YTickRange_Sh,'yTickLabel',YTickLabel_Sh,...
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
clearvars -except Met_Year
end
