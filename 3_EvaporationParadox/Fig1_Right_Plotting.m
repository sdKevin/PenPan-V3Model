function Fig1_Right_Plotting(Epan_Year,YlimRange_Epan,YlimRange_Epan_R,YlimRange_Epan_A)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Epan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];

%% Plotting Epan Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % E_pan m/s to mm/year
    Epan_Historical = Epan_Year(1).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    Epan_ssp(:,:,i_ssp-1) = Epan_Year(i_ssp).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    Epan_ssp(:,:,i_ssp-1) = Epan_ssp(:,:,i_ssp-1) - repmat(mean(Epan_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_ssp(:,end-30:end-1,1)')',mean(Epan_ssp(:,end-30:end-1,2)')',...
    mean(Epan_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Epan Anomaly (mm year^-^1)');
%% Plotting Epan Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
     % E_pan m/s to mm/year
    Epan_Historical = Epan_Year(1).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    Epan_Historical(16,:) = []; % HadGEM3-GC31-LL
    Epan_ssp370 = Epan_Year(i_ssp).Epan_Year.E_pan .* 365.*24.*3600.*1000;
    Epan_ssp370 = Epan_ssp370 - repmat(mean(Epan_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_ssp370(:,end-30:end-1)')',mean(Epan_ssp370(:,end-30:end-1)')',...
    mean(Epan_ssp370(:,end-30:end-1)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Epan_Year YlimRange_Epan YlimRange_Epan_R YlimRange_Epan_A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Epan_R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];

%% Plotting Epan_R Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Epan_R m/s to mm/year
    Epan_R_Historical = Epan_Year(1).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    Epan_R_ssp(:,:,i_ssp-1) = Epan_Year(i_ssp).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    Epan_R_ssp(:,:,i_ssp-1) = Epan_R_ssp(:,:,i_ssp-1) - repmat(mean(Epan_R_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_R_ssp(:,end-30:end-1,1)')',mean(Epan_R_ssp(:,end-30:end-1,2)')',...
    mean(Epan_R_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan_R,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Epan_R Anomaly (mm year^-^1)');
%% Plotting Epan_R Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Epan_R m/s to mm/year
    Epan_R_Historical = Epan_Year(1).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    Epan_R_Historical(16,:) = []; % HadGEM3-GC31-LL
    Epan_R_ssp370 = Epan_Year(i_ssp).Epan_Year.E_pan_R .* 365.*24.*3600.*1000;
    Epan_R_ssp370 = Epan_R_ssp370 - repmat(mean(Epan_R_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_R_ssp370(:,end-30:end-1)')',mean(Epan_R_ssp370(:,end-30:end-1)')',...
    mean(Epan_R_ssp370(:,end-30:end-1)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan_R,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Epan_Year YlimRange_Epan YlimRange_Epan_R YlimRange_Epan_A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Epan_A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,6,[1:3])
%% Setting Color and axis property
RGB_ssp_Shade = [222,235,247; 198,234,251; 161,196,218; 194,196,226];
RGB_ssp_Line = [133,184,227; 0,173,238; 50,128,185; 57,83,164];

%% Plotting Epan_A Under ssp126 , 245 & 585
for i_ssp = [5 , 3 , 2]
    % Epan_A m/s to mm/year
    Epan_A_Historical = Epan_Year(1).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    Epan_A_ssp(:,:,i_ssp-1) = Epan_Year(i_ssp).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    Epan_A_ssp(:,:,i_ssp-1) = Epan_A_ssp(:,:,i_ssp-1) - repmat(mean(Epan_A_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_A_ssp(:,end-30:end-1,1)')',mean(Epan_A_ssp(:,end-30:end-1,2)')',...
    mean(Epan_A_ssp(:,end-30:end-1,4)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan_A,'FontSize',24,'FontName','Arial','LineWidth',2.5);
ylabel('Epan_A Anomaly (mm year^-^1)');
%% Plotting Epan_A Under ssp370
subplot(1,6,[4:6])
for i_ssp = 4
    % Epan_A m/s to mm/year
    Epan_A_Historical = Epan_Year(1).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    Epan_A_Historical(16,:) = []; % HadGEM3-GC31-LL
    Epan_A_ssp370 = Epan_Year(i_ssp).Epan_Year.E_pan_A .* 365.*24.*3600.*1000;
    Epan_A_ssp370 = Epan_A_ssp370 - repmat(mean(Epan_A_Historical(:,99:165),2), 1 ,86); % Change to Anomaly
end

H = notBoxPlot([mean(Epan_A_ssp370(:,end-30:end-1)')',mean(Epan_A_ssp370(:,end-30:end-1)')',...
    mean(Epan_A_ssp370(:,end-30:end-1)')'],'jitter',0.3);
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

set(gca,'ylim',YlimRange_Epan_A,'FontSize',24,'FontName','Arial','LineWidth',2.5);
clearvars -except Epan_Year YlimRange_Epan YlimRange_Epan_R YlimRange_Epan_A
end