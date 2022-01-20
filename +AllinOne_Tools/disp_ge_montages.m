function [Anatomic_Fig,Mask_Fig,VentMontage,GasMontage,DissolvedMontage,RBCMontage,BarrierMontage,VentBinMontage,DissolvedBinMontage,RBCBinMontage,BarrierBinMontage,RBCBarBinMontage] = disp_ge_montages(H1_Image_Dis,Proton_Mask,ScaledVentImage,LoRes_Gas_Image,Dis_Image,Bar_Image,RBC_Image,RBC2BarIm,VentBinMap,DissolvedBinMap,BarrierBinMap,RBCBinMap,RBCBarrierBinMap,SNRS);

ProtonMax = prctile(abs(H1_Image_Dis(:)),99.99);

VentSNR = SNRS.VentSNR;
GasSNR = SNRS.GasSNR;
DissolvedSNR = SNRS.DissolvedSNR;
BarrierSNR = SNRS.BarrierSNR;
RBCSNR = SNRS.RBCSNR;
RBC2Bar = SNRS.RBC2Bar;

SixBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 0 0.57 0.71; 0 0 1]; %Used for Vent and RBC
EightBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255]; %Used for barrier
SixBinRBCBarMap = [ 197/255 27/255 125/255; 225/255 129/255 162/255; 0.4 0.7 0.4; 0 1 0; 1 0.7143 0; 1 0 0]; %Used for RBC/Barrier ratio

%% Display all the images - Do this in a better way than just using montage
Vent_Dis_RBC_Label = {'Defect','Low','Healthy','Healthy','High','High'};
Bar_Label = {'Defect','Low','Healthy','Healthy','Elevated','Elevated','High','High'};
RBCBar_Label = {'Bar Dom','Bar Dom','Healthy','Healthy','RBC Dom','RBC Dom'};

[~,firstslice,lastslice] = AllinOne_Tools.getimcenter(Proton_Mask);
%Make a bunch of tiled images for making nice figures
try
    Anat_tiled1 = AllinOne_Tools.tile_image(H1_Image_Dis(:,:,(firstslice-2):(lastslice+2)),3,'nRows',3);
    Mask_tiled = AllinOne_Tools.tile_image(Proton_Mask(:,:,(firstslice-2):(lastslice+2)),3,'nRows',3);
    Vent_tiled = AllinOne_Tools.tile_image(ScaledVentImage(:,:,(firstslice-2):(lastslice+2)),3,'nRows',3);
    Gas_tiled = AllinOne_Tools.tile_image(LoRes_Gas_Image(:,:,(firstslice-2):(lastslice+2)),3,'nRows',3);
    Dissolved_tiled = AllinOne_Tools.tile_image(Dis_Image(:,:,(firstslice-2):(lastslice+2)),3,'nRows',3);
    Barrier_tiled = AllinOne_Tools.tile_image(abs(Bar_Image(:,:,(firstslice-2):(lastslice+2))),3,'nRows',3);
    RBC_tiled = AllinOne_Tools.tile_image(abs(RBC_Image(:,:,(firstslice-2):(lastslice+2))),3,'nRows',3);
    Anat_tiled = AllinOne_Tools.tile_image(H1_Image_Dis(:,:,(firstslice-2):(lastslice+2)),3);
    RBCBar_tiled = AllinOne_Tools.tile_image(RBC2BarIm(:,:,(firstslice-2):(lastslice+2)),3);
    VentBin_tiled = AllinOne_Tools.tile_image(VentBinMap(:,:,(firstslice-2):(lastslice+2)),3);
    DissolvedBin_tiled = AllinOne_Tools.tile_image(DissolvedBinMap(:,:,(firstslice-2):(lastslice+2)),3);
    BarrierBin_tiled = AllinOne_Tools.tile_image(BarrierBinMap(:,:,(firstslice-2):(lastslice+2)),3);
    RBCBin_tiled = AllinOne_Tools.tile_image(RBCBinMap(:,:,(firstslice-2):(lastslice+2)),3);
    RBCBarBin_tiled = AllinOne_Tools.tile_image(RBCBarrierBinMap(:,:,(firstslice-2):(lastslice+2)),3);
catch
    Anat_tiled1 = AllinOne_Tools.tile_image(H1_Image_Dis,3,'nRows',3);
    Mask_tiled = AllinOne_Tools.tile_image(Mask,3,'nRows',3);
    Vent_tiled = AllinOne_Tools.tile_image(ScaledVentImage,3,'nRows',3);
    Gas_tiled = AllinOne_Tools.tile_image(LoRes_Gas_Image,3,'nRows',3);
    Dissolved_tiled = AllinOne_Tools.tile_image(Dis_Image,3,'nRows',3);
    Barrier_tiled = AllinOne_Tools.tile_image(Bar_Image,3,'nRows',3);
    RBC_tiled = AllinOne_Tools.tile_image(RBC_Image,3,'nRows',3);
    Anat_tiled = AllinOne_Tools.tile_image(H1_Image_Dis,3);
    RBCBar_tiled = AllinOne_Tools.tile_image(RBC2BarIm,3);
    VentBin_tiled = AllinOne_Tools.tile_image(VentBinMap,3);
    DissolvedBin_tiled = AllinOne_Tools.tile_image(DissolvedBinMap,3);
    BarrierBin_tiled = AllinOne_Tools.tile_image(BarrierBinMap,3);
    RBCBin_tiled = AllinOne_Tools.tile_image(RBCBinMap,3);
    RBCBarBin_tiled = AllinOne_Tools.tile_image(RBCBarrierBinMap,3);
end

Anatomic_Fig = figure('Name','All Anatomic','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(Anatomic_Fig,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(Anatomic_Fig,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Anat_tiled1,[0.9,1.1],[0,0.99*ProtonMax],gray,0,gca);
axis off
colormap(gray);
title('Anatomic Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
set(Anatomic_Fig,'WindowState','minimized');

Mask_Map = [1 0 0];
Mask_Fig = figure('Name','All Anatomic with Mask','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(Mask_Fig,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(Mask_Fig,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Mask_tiled,[0.9,1.1],[0,0.99*ProtonMax],Mask_Map,0.25,gca);
colormap(gca,Mask_Map)
title('Anatomic Image Masked','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
set(Mask_Fig,'WindowState','minimized');

Scaled_Vent_tile = Vent_tiled/(max(Vent_tiled(Mask_tiled==1)));
%Plot Ventilation with Mask outlines
VentMontage = figure('Name','All Ventilation with Mask Outline','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(VentMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(VentMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Scaled_Vent_tile,[0,1],[0,0.99*ProtonMax],gray,1,gca);
%imagesc(Vent_tiled);
%axis off
colormap(gray)
hold on
B = bwboundaries(Mask_tiled);
for j = 1:length(B)
    if length(B{j})>2
        plot(B{j}(:,2),B{j}(:,1),'r')
    end
end
title('Ventilation Image with Mask Boundaries','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(VentMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(VentSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(VentMontage,'WindowState','minimized');

GasMontage = figure('Name','Low Resolution Gas','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(GasMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(GasMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Gas_tiled,[0 max(abs(LoRes_Gas_Image(:)))],[0 max(abs(LoRes_Gas_Image(:)))],gray,1,gca);
axis off
colormap(gray);
title('Low Resolution Gas Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(GasMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(GasSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(GasMontage,'WindowState','minimized');

DissolvedMontage = figure('Name','Dissolved','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(DissolvedMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(DissolvedMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Dissolved_tiled,[0 max(abs(Dissolved_tiled(:)))],[0 max(abs(Dissolved_tiled(:)))],gray,1,gca);
axis off
colormap(gray);
title('Dissolved Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(DissolvedMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(DissolvedSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(DissolvedMontage,'WindowState','minimized');

RBCMontage = figure('Name','RBC','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(RBCMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(RBCMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,RBC_tiled,[0 max(abs(RBC_tiled(:)))],[0 max(abs(RBC_tiled(:)))],gray,1,gca);
axis off
colormap(gray);
title('RBC Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(RBCMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(RBCSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(RBCMontage,'WindowState','minimized');

BarrierMontage = figure('Name','Barrier','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
set(BarrierMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
%set(BarrierMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled1,Barrier_tiled,[0 max(abs(Barrier_tiled(:)))],[0 max(abs(Barrier_tiled(:)))],gray,1,gca);
axis off
colormap(gray);
title('Barrier Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(BarrierMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(BarrierSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(BarrierMontage,'WindowState','minimized');

% Now Binned Montages
VentBinMontage = figure('Name','Binned Ventilation','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
%set(VentBinMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
set(VentBinMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled,VentBin_tiled,[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
colormap(gca,SixBinMap)
cbar = colorbar(gca','Location','southoutside','Ticks',[]);
pos = cbar.Position;
cbar.Position = [pos(1),0,pos(3),pos(4)];
try
    AllinOne_Tools.binning_colorbar(cbar,6,Vent_Dis_RBC_Label);
catch
end
title('Binned Ventilation Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(VentBinMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(VentSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(VentBinMontage,'WindowState','minimized');

DissolvedBinMontage = figure('Name','Dissolved','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
%set(DissolvedBinMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
set(DissolvedBinMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled,DissolvedBin_tiled,[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
colormap(gca,SixBinMap)
cbar = colorbar(gca','Location','southoutside','Ticks',[]);
pos = cbar.Position;
cbar.Position = [pos(1),0,pos(3),pos(4)];
try
    AllinOne_Tools.binning_colorbar(cbar,6,Vent_Dis_RBC_Label);
catch
end
axis off
title('Binned Dissolved Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(DissolvedBinMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(DissolvedSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(DissolvedBinMontage,'WindowState','minimized');

RBCBinMontage = figure('Name','RBC','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
%set(RBCBinMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
set(RBCBinMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled,RBCBin_tiled,[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
axis off
colormap(gca,SixBinMap)
cbar = colorbar(gca','Location','southoutside','Ticks',[]);
pos = cbar.Position;
cbar.Position = [pos(1),0,pos(3),pos(4)];
try
    AllinOne_Tools.binning_colorbar(cbar,6,Vent_Dis_RBC_Label);
catch
end
title('Binned RBC Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(RBCBinMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(RBCSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(RBCBinMontage,'WindowState','minimized');

BarrierBinMontage = figure('Name','Barrier','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
%set(BarrierBinMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
set(BarrierBinMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled,BarrierBin_tiled,[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
axis off
colormap(gca,EightBinMap)
cbar = colorbar(gca','Location','southoutside','Ticks',[]);
pos = cbar.Position;
cbar.Position = [pos(1),0,pos(3),pos(4)];
try
    AllinOne_Tools.binning_colorbar(cbar,8,Bar_Label);
catch
end
title('Binned Barrier Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(BarrierBinMontage,'textbox',[0.8 0.08 0.2 0.05],'Color',[1 1 1],'String',['SNR = ' num2str(BarrierSNR,'%.1f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(BarrierBinMontage,'WindowState','minimized');

RBCBarBinMontage = figure('Name','RBC to Barrier','units','normalized','outerposition',[.2 .2 1 4/3]);%set(ClinFig,'WindowState','minimized');
%set(RBCBarBinMontage,'color','white','Units','inches','Position',[1 1 10 3.3])
set(RBCBarBinMontage,'color','white','Units','inches','Position',[1 1 8 7.2])
axes('Units', 'normalized', 'Position', [0 0 1 1])
[~,~] = AllinOne_Tools.imoverlay(Anat_tiled,RBCBarBin_tiled,[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
axis off
colormap(gca,SixBinRBCBarMap)
cbar = colorbar(gca','Location','southoutside','Ticks',[]);
pos = cbar.Position;
cbar.Position = [pos(1),0,pos(3),pos(4)];
try
    AllinOne_Tools.binning_colorbar(cbar,6,RBCBar_Label);
catch
end
title('Binned RBC/Barrier Image','FontSize',16)
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)-.01])
annotation(RBCBarBinMontage,'textbox',[0.7 0.08 0.2 0.05],'Color',[1 1 1],'String',['RBC/Barrier = ' num2str(RBC2Bar,'%.2f')],'FontSize',14,'FontName','Arial','FitBoxToText','on','BackgroundColor',[0 0 0],'VerticalAlignment','middle','HorizontalAlignment','center');
set(RBCBarBinMontage,'WindowState','minimized');

