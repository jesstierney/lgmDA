% Written by M. Osman, Jan2020 (mattosman@arizona.edu) in Matlab2018b
% This function plots the integral timescale of TAS, SOS and TOS for each
% iCESM LGM "slice" (PI, 3ky BP, 6ky BP, 9ky BP, 12ky BP, 18kyBP, 21kyBP)
% 
% This function requires:
% 1. A folder named ''ModelData/'', housing the following files:
%   a. d18Osw.mat		
%   b. toGA.mat
%   c. sos.mat			
%   d. tos.mat
%   e. tas.mat	
% 2. the function integtimsc.m
clear all
home_fold = cd;

cd ModelData/
	load d18Osw.mat		
	load toGA.mat
	load sos.mat			
	load tos.mat
    latPOP = double(lat); lonPOP = double(lon);
    load tas.mat	
    lat = double(lat);  lon = double(lon);  
cd ../

% compute for each slice, for tas, tos, and sos
tas_intTS = nan(size(tas,1),size(tas,2),size(tas,4));
tos_intTS = nan(size(tos,1),size(tos,2),size(tos,4));
sos_intTS = nan(size(sos,1),size(sos,2),size(sos,4));
% atmos
for i = 1:size(tas,4)
    for j = 1:size(tas,1)
        for k = 1:size(tas,2)
            tas_intTS(j,k,i) = integtimsc(squeeze(tas(j,k,:,i)))*10;
        end
    end
end

% ocean
for i = 1:size(tos,4)
    for j = 1:size(tos,1)
        for k = 1:size(tos,2)
            tos_intTS(j,k,i) = integtimsc(squeeze(tos(j,k,:,i)))*10;
            sos_intTS(j,k,i) = integtimsc(squeeze(sos(j,k,:,i)))*10;
        end
    end
end

%% TAS Figure;

h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.55*Pix_SS(3),.70*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,101); 
    CT2 = cbrewer('seq','Blues' ,101);
warning('on','all')
cd(home_fold)

load coastlines.mat

% find global max min vals
for i = 1:size(tas_intTS,3)
    currMap = squeeze(tas_intTS(:,:,i)); 
    if i == 1
        globMax = nanmax(currMap(:));
        globMin = nanmin(currMap(:));
    else
        if nanmax(currMap(:)) > globMax
            globMax = nanmax(currMap(:));
        elseif nanmin(currMap(:)) < globMin
            globMin = nanmin(currMap(:));
        end
    end
end    

% PI
currMap = squeeze(tas_intTS(:,:,1));
ax1 = axes('Position',[0.20 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');
       
% 3ky
currMap = squeeze(tas_intTS(:,:,2));
ax1 = axes('Position',[0.20 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');
    
% 6ky
currMap = squeeze(tas_intTS(:,:,3));
ax1 = axes('Position',[0.20 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');    
    
% 9ky
currMap = squeeze(tas_intTS(:,:,4));
ax1 = axes('Position',[0.20 0.0 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');   
    
% 12ky
currMap = squeeze(tas_intTS(:,:,5));
ax1 = axes('Position',[0.50 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');   
    
% 18ky
currMap = squeeze(tas_intTS(:,:,6));
ax1 = axes('Position',[0.50 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');   
    
% 21ky
currMap = squeeze(tas_intTS(:,:,7));
ax1 = axes('Position',[0.50 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(lat,lon,currMap'); shading flat;
    colormap(ax1,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    plotm(coastlat,coastlon,'-k');     
    
    
ax5 = axes('position', [0.05 0.05 0.10 0.10]); hold on; box off;
set(gca, 'visible', 'off')
text(0,0,'PI','Fontsize',11,'Fontweight','normal');
text(0,0,'3-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'6-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'9-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'12-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'18-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'21-ky BP','Fontsize',11,'Fontweight','normal');

%% TOS Figure

h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.55*Pix_SS(3),.70*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,101); 
    CT2 = cbrewer('seq','Blues' ,101);
warning('on','all')
cd(home_fold)

land = shaperead('landareas.shp','UseGeoCoords',true);

% find global max min vals
for i = 1:size(tos_intTS,3)
    currMap = squeeze(tos_intTS(:,:,i)); 
    if i == 1
        globMax = nanmax(currMap(:));
        globMin = nanmin(currMap(:));
    else
        if nanmax(currMap(:)) > globMax
            globMax = nanmax(currMap(:));
        elseif nanmin(currMap(:)) < globMin
            globMin = nanmin(currMap(:));
        end
    end
end    

% PI
currMap = squeeze(tos_intTS(:,:,1));
ax1 = axes('Position',[0.20 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
       
% 3ky
currMap = squeeze(tos_intTS(:,:,2));
ax1 = axes('Position',[0.20 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 6ky
currMap = squeeze(tos_intTS(:,:,3));
ax1 = axes('Position',[0.20 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 9ky
currMap = squeeze(tos_intTS(:,:,4));
ax1 = axes('Position',[0.20 0.0 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 12ky
currMap = squeeze(tos_intTS(:,:,5));
ax1 = axes('Position',[0.50 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 18ky
currMap = squeeze(tos_intTS(:,:,6));
ax1 = axes('Position',[0.50 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 21ky
currMap = squeeze(tos_intTS(:,:,7));
ax1 = axes('Position',[0.50 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
ax5 = axes('position', [0.05 0.05 0.10 0.10]); hold on; box off;
set(gca, 'visible', 'off')
text(0,0,'PI','Fontsize',11,'Fontweight','normal');
text(0,0,'3-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'6-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'9-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'12-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'18-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'21-ky BP','Fontsize',11,'Fontweight','normal');

%% SOS Figure

h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.55*Pix_SS(3),.70*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,101); 
    CT2 = cbrewer('seq','Blues' ,101);
    CT3 = cbrewer('seq','Purples' ,101);
    % CT3 = cbrewer('div','RdBu' ,101);
warning('on','all')
cd(home_fold)

land = shaperead('landareas.shp','UseGeoCoords',true);

% find global max min vals
for i = 1:size(sos_intTS,3)
    currMap = squeeze(sos_intTS(:,:,i)); 
    if i == 1
        globMax = nanmax(currMap(:));
        globMin = nanmin(currMap(:));
    else
        if nanmax(currMap(:)) > globMax
            globMax = nanmax(currMap(:));
        elseif nanmin(currMap(:)) < globMin
            globMin = nanmin(currMap(:));
        end
    end
end    

% PI
currMap = squeeze(sos_intTS(:,:,1));
ax1 = axes('Position',[0.20 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
       
% 3ky
currMap = squeeze(sos_intTS(:,:,2));
ax1 = axes('Position',[0.20 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 6ky
currMap = squeeze(sos_intTS(:,:,3));
ax1 = axes('Position',[0.20 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 9ky
currMap = squeeze(sos_intTS(:,:,4));
ax1 = axes('Position',[0.20 0.0 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 12ky
currMap = squeeze(sos_intTS(:,:,5));
ax1 = axes('Position',[0.50 0.675 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 18ky
currMap = squeeze(sos_intTS(:,:,6));
ax1 = axes('Position',[0.50 0.45 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
% 21ky
currMap = squeeze(sos_intTS(:,:,7));
ax1 = axes('Position',[0.50 0.225 0.325 0.325]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12); % mlabel(0); 
    pcolorm(latPOP,lonPOP,currMap); shading flat;
    colormap(ax1,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'Memory (years)'}]; set(c,'Fontsize',11);
        caxis([globMin globMax]);
        % caxis([0 200]);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    
ax5 = axes('position', [0.05 0.05 0.10 0.10]); hold on; box off;
set(gca, 'visible', 'off')
text(0,0,'PI','Fontsize',11,'Fontweight','normal');
text(0,0,'3-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'6-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'9-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'12-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'18-ky BP','Fontsize',11,'Fontweight','normal');
text(0,0,'21-ky BP','Fontsize',11,'Fontweight','normal');
