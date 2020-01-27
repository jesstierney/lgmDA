%% Plot sst proxies by type and location, and sampling frequency
% written by m. osman Nov 2019 (mattosman@arizona.edu)
clear
home_fold = cd;
BinSize = 250; % sampling frequency in years, for plotting

% load in the proxies;
cd sst_proxy_repository
    load uk37_proxies_sorted_05-Nov-2019.mat % uk37
    load tex86_proxies_sorted_05-Nov-2019.mat % tex86
    load d18o_proxies_sorted_05-Nov-2019.mat % d18O
    load mgca_proxies_sorted_05-Nov-2019.mat % mgca
cd ../

% Define figure parameters 
h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.55*Pix_SS(3),.70*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,5); 
    CT2 = cbrewer('seq','Blues' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);    
    CT5 = cbrewer('div','Spectral' ,length(1850:10:2010));    
    cd '/Users/matthewosman/Documents/MATLAB/colornames';
    [~,gold] = colornames('CSS','Goldenrod');
warning('on','all')
cd(home_fold)

ax1 = axes('Position',[0.15 0.30 0.70 0.80]); 
    set(gca,'Visible','off'); 
    ax = axesm ('robinson','MLabelLocation',60, 'Frame', 'on', 'Grid', 'off');
        gridm('mlinelocation',[],'plinelocation',90,'mlinelocation',60); 
        ax.Clipping = 'off'; set(ax,'Fontsize',12)
    land = shaperead('landareas.shp','UseGeoCoords',true);
    geoshow(land,'FaceColor', [0.4 0.4 0.4],'Linewidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
    % v1
    for i = 1:length(d18o.data)
        a1 = scatterm(d18o.lats(i),d18o.lons(i), 150, CT2(end-1,:),'Filled'); a1.Children.Marker = 'o'; % 'square';         
            a1.Children.MarkerEdgeColor = 'none'; a1.Children.MarkerFaceAlpha = 0.90; % 0.45;
    end    
    for i = 1:length(uk37.data)
        a2 = scatterm(uk37.lats(i),uk37.lons(i), 120, CT1(end-1,:),'Filled'); a2.Children.Marker = 'o'; % '>';  
            a2.Children.MarkerEdgeColor = 'none'; a2.Children.MarkerFaceAlpha = 0.90; % 0.45;
    end
    for i = 1:length(mgca.data)
        a3 = scatterm(mgca.lats(i),mgca.lons(i), 90, gold,'Filled'); a3.Children.Marker = 'o'; % 'o';         
            a3.Children.MarkerEdgeColor = 'none'; a3.Children.MarkerFaceAlpha = 0.90; % 0.45;
    end    
    for i = 1:length(tex86.data)
        a4 = scatterm(tex86.lats(i),tex86.lons(i), 60, CT3(end-1,:),'Filled'); a4.Children.Marker = 'o'; % 'diamond';
            a4.Children.MarkerEdgeColor = 'none'; a4.Children.MarkerFaceAlpha = 0.90; % 0.45;
    end    
    legend([a1,a2,a3,a4],'\delta_{18}O_{foram.}', 'U^{K''}_{37}', 'MgCa', 'TEX_{86}', 'Orientation','Vertical','Fontsize',11,'Box','off')
%     % v2
%     for i = 1:length(d18o.data)
%         a1 = scatterm(d18o.lats(i),d18o.lons(i), 50, CT2(end-1,:));
%             a1.Children.LineWidth = 0.5; a1.Children.Marker = 'square';         
%     end    
%     for i = 1:length(uk37.data)
%         a2 = scatterm(uk37.lats(i),uk37.lons(i), 50, CT1(end-1,:));
%             a2.Children.LineWidth = 0.5; a2.Children.Marker = '>';  
%     end
%     for i = 1:length(mgca.data)
%         a3 = scatterm(mgca.lats(i),mgca.lons(i), 50, gold);
%             a3.Children.LineWidth = 0.5; a3.Children.Marker = 'o';         
%     end    
%     for i = 1:length(tex86.data)
%         a4 = scatterm(tex86.lats(i),tex86.lons(i), 50, CT3(end-1,:));
%             a4.Children.LineWidth = 0.5; a4.Children.Marker = 'diamond';
%     end    
%     legend([a1,a2,a3,a4],'\delta_{18}O_{foram.}', 'U^{K''}_{37}', 'MgCa', 'TEX_{86}', 'Orientation','Vertical','Fontsize',11,'Box','off')

    
    % subplot2: number of datapoints per bin size;
    ageBin = [0:BinSize:25000]';
    % step1 = concatenate all datapoints for each proxy type:
    d18O_allAges  = [];
    mgca_allAges  = [];
    uk37_allAges  = [];
    tex86_allAges = [];
    for i = 1:length(d18o.age); d18O_allAges = [d18O_allAges; d18o.age{i}]; end
    for i = 1:length(mgca.age); mgca_allAges = [mgca_allAges; mgca.age{i}]; end
    for i = 1:length(uk37.age); uk37_allAges = [uk37_allAges; uk37.age{i}]; end
    for i = 1:length(tex86.age); tex86_allAges = [tex86_allAges; tex86.age{i}]; end
    % now, count the number of data values whose age falls within each bin
    d18O_binNum = zeros(length(ageBin)-1,1);
    mgca_binNum = zeros(length(ageBin)-1,1);
    uk37_binNum = zeros(length(ageBin)-1,1);
    tex86_binNum = zeros(length(ageBin)-1,1);
    for i = 1:length(ageBin)-1
        indexer = d18O_allAges >= ageBin(i) & d18O_allAges < ageBin(i+1);
        d18O_binNum(i) = nansum(indexer); clearvars indexer;
    end
    for i = 1:length(ageBin)-1
        indexer = mgca_allAges >= ageBin(i) & mgca_allAges < ageBin(i+1);
        mgca_binNum(i) = nansum(indexer); clearvars indexer;
    end    
    for i = 1:length(ageBin)-1
        indexer = uk37_allAges >= ageBin(i) & uk37_allAges < ageBin(i+1);
        uk37_binNum(i) = nansum(indexer); clearvars indexer;
    end    
    for i = 1:length(ageBin)-1
        indexer = tex86_allAges >= ageBin(i) & tex86_allAges < ageBin(i+1);
        tex86_binNum(i) = nansum(indexer); clearvars indexer;
    end   
    ageBin = (ageBin(1:end-1) + ageBin(2:end))./2; % adjust so mid point for plotting
    combined = [d18O_binNum, mgca_binNum, uk37_binNum, tex86_binNum];
    % plot
    ax2 = axes('Position',[0.15 0.10 0.70 0.30]); 
    b = bar(ageBin,combined,0.95,'stacked','EdgeColor','none');
        b(1).FaceColor = CT2(end-1,:);
        b(2).FaceColor = CT1(end-1,:);
        b(3).FaceColor = gold;
        b(4).FaceColor = CT3(end-1,:);
    xlabel('Age BP'); ylabel('Proxy samples bin^{-1}')
    set(gca,'color','none','yaxislocation','left','linewidth',1.5,'Fontsize',12); grid on; box off;
    set(gca,'Xtick',[0,5000,10000,15000,20000,25000],'Xticklabel',...
        [{'0'},{'5000'},{'10000'},{'15000'},{'20000'},{'25000'}])