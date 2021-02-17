function Fig2_Right_Plotting(Met_Year)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Sg = [-6,6];

%% Plotting Sg Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Sg W/m2
    Sg_Historical = Met_Year(1).Met_Year.Sg;
    Sg_ssp(:,:,i_ssp-1) = Met_Year(i_ssp).Met_Year.Sg;
    Sg_ssp(:,:,i_ssp-1) = Sg_ssp(:,:,i_ssp-1) - repmat(mean(Sg_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Sg_ssp(:,end-30:end-1,1)')',mean(Sg_ssp(:,end-30:end-1,2)')',...
    mean(Sg_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp126';'ssp245';'ssp585']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(1,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(1,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(2,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(2,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(4,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(4,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Sg,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Sg Anomaly (W m^-^2)');
%% Plotting Sg Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Sg W/m2
    Sg_Historical = Met_Year(1).Met_Year.Sg;
    Sg_Historical(16,:) = []; % HadGEM3-GC31-LL
    Sg_ssp370 = Met_Year(i_ssp).Met_Year.Sg;
    Sg_ssp370 = Sg_ssp370 - repmat(mean(Sg_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Sg_ssp370(:,end-30:end-1)')',mean(Sg_ssp370(:,end-30:end-1)')',...
    mean(Sg_ssp370(:,end-30:end-1)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp370','','']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Sg,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Li %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Li = [-10,50];

%% Plotting Li Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Li W/m2
    Li_Historical = Met_Year(1).Met_Year.Li;
    Li_ssp(:,:,i_ssp-1) = Met_Year(i_ssp).Met_Year.Li;
    Li_ssp(:,:,i_ssp-1) = Li_ssp(:,:,i_ssp-1) - repmat(mean(Li_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Li_ssp(:,end-30:end-1,1)')',mean(Li_ssp(:,end-30:end-1,2)')',...
    mean(Li_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp126';'ssp245';'ssp585']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(1,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(1,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(2,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(2,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(4,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(4,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Li,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Li Anomaly (W m^-^2)');
%% Plotting Li Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Li W/m2
    Li_Historical = Met_Year(1).Met_Year.Li;
    Li_Historical(16,:) = []; % HadGEM3-GC31-LL
    Li_ssp370 = Met_Year(i_ssp).Met_Year.Li;
    Li_ssp370 = Li_ssp370 - repmat(mean(Li_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Li_ssp370(:,end-30:end-1)')',mean(Li_ssp370(:,end-30:end-1)')',...
    mean(Li_ssp370(:,end-30:end-1)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp370','','']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Li,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% U10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_U10 = [-0.08,0.06];

%% Plotting U10 Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % U10 m/s
    U10_Historical = Met_Year(1).Met_Year.U10;
    U10_ssp(:,:,i_ssp-1) = Met_Year(i_ssp).Met_Year.U10;
    U10_ssp(:,:,i_ssp-1) = U10_ssp(:,:,i_ssp-1) - repmat(mean(U10_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(U10_ssp(:,end-30:end-1,1)')',mean(U10_ssp(:,end-30:end-1,2)')',...
    mean(U10_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp126';'ssp245';'ssp585']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(1,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(1,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(2,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(2,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(4,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(4,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_U10,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('U10 Anomaly (m s^-^1)');
%% Plotting U10 Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % U10 m/s
    U10_Historical = Met_Year(1).Met_Year.U10;
    U10_Historical(16,:) = []; % HadGEM3-GC31-LL
    U10_ssp370 = Met_Year(i_ssp).Met_Year.U10;
    U10_ssp370 = U10_ssp370 - repmat(mean(U10_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(U10_ssp370(:,end-30:end-1)')',mean(U10_ssp370(:,end-30:end-1)')',...
    mean(U10_ssp370(:,end-30:end-1)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp370','','']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_U10,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Ta = [-2,8];

%% Plotting Ta Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Ta K
    Ta_Historical = Met_Year(1).Met_Year.Ta;
    Ta_ssp(:,:,i_ssp-1) = Met_Year(i_ssp).Met_Year.Ta;
    Ta_ssp(:,:,i_ssp-1) = Ta_ssp(:,:,i_ssp-1) - repmat(mean(Ta_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Ta_ssp(:,end-30:end-1,1)')',mean(Ta_ssp(:,end-30:end-1,2)')',...
    mean(Ta_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp126';'ssp245';'ssp585']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(1,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(1,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(2,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(2,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(4,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(4,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Ta,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Ta Anomaly (K)');
%% Plotting Ta Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Ta K
    Ta_Historical = Met_Year(1).Met_Year.Ta;
    Ta_Historical(16,:) = []; % HadGEM3-GC31-LL
    Ta_ssp370 = Met_Year(i_ssp).Met_Year.Ta;
    Ta_ssp370 = Ta_ssp370 - repmat(mean(Ta_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Ta_ssp370(:,end-30:end-1)')',mean(Ta_ssp370(:,end-30:end-1)')',...
    mean(Ta_ssp370(:,end-30:end-1)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp370','','']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Ta,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Met_Year

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];
YlimRange_Sh = [-0.5,3.5]./1000;

%% Plotting Sh Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Sh (-)
    Sh_Historical = Met_Year(1).Met_Year.Sh;
    Sh_ssp(:,:,i_ssp-1) = Met_Year(i_ssp).Met_Year.Sh;
    Sh_ssp(:,:,i_ssp-1) = Sh_ssp(:,:,i_ssp-1) - repmat(mean(Sh_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Sh_ssp(:,end-30:end-1,1)')',mean(Sh_ssp(:,end-30:end-1,2)')',...
    mean(Sh_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp126';'ssp245';'ssp585']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(1,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(1,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(2,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(2,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'Color',1.0-1*(1.0-RGB_ssp_Line(4,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(4,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Sh,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Sh Anomaly (¡Á10^-^3)');
%% Plotting Sh Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Sh (-)
    Sh_Historical = Met_Year(1).Met_Year.Sh;
    Sh_Historical(16,:) = []; % HadGEM3-GC31-LL
    Sh_ssp370 = Met_Year(i_ssp).Met_Year.Sh;
    Sh_ssp370 = Sh_ssp370 - repmat(mean(Sh_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Sh_ssp370(:,end-30:end-1)')',mean(Sh_ssp370(:,end-30:end-1)')',...
    mean(Sh_ssp370(:,end-30:end-1)')'],'jitter',0.3);
set(gca,'color',[ 255 255 255]./255,'Box','On','xTickLabel',['ssp370','','']);
set([H.data],'MarkerSize',6,'markerFaceColor','none','markerEdgeColor', 'none');

set([H(1).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(1).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(1).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(2).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(2).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(2).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set([H(3).sdPtch],'EdgeColor','none','FaceAlpha',0); %delete 1std
set([H(3).mu],'LineStyle',':','Color',1.0-1*(1.0-RGB_ssp_Line(3,:)./255),'EraseMode','xor');
set([H(3).semPtch],'FaceColor',RGB_ssp_Shade(3,:)./255,'FaceAlpha',0.9,'EdgeColor','none');

set(gca,'ylim',YlimRange_Sh,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Met_Year
end