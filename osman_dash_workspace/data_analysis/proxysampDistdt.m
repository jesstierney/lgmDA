function [d18o_dist,mgca_dist,uk37_dist,tex86_dist,ageBin] = proxysampDistdt(stepSize,windSize,plotOpt)

% Written by M. Osman, Jan2020 (mattosman@arizona.edu) in Matlab2018b
% This function determines the distribution of sampling frequencies (i.e.,
% resolutions) of the LGM - present DA proxy compilation (see Tierney et al., 
% in review).  Insight into the resolution (including potential resolution
% time-dependency) is an important consideration for determining the
% binning size of proxies prior to assimilating. 
% This function follows from proxysampDist.m but analyzes change in sample
% res distributions through time
% 
% This function requires:
% 1. A folder named ''sst_proxy_repository/'', housing the following files:
%   a. d18o_proxies_sorted_21-Jan-2020.mat
%   b. mgca_proxies_sorted_21-Jan-2020.mat
%   c. tex86_proxies_sorted_21-Jan-2020.mat
%   d. uk37_proxies_sorted_21-Jan-2020.mat
% NOTE these files can be created via extractProx.m" 
home_fold = cd;

% check stepSize input
if ~isnumeric(stepSize)
    warning('Input for ''binSize'' must be numeric. Try again.'); return;
elseif stepSize <= 0 || stepSize > 1000
    warning('Input for ''binSize'' must be >= 0 or < 1000. Try again.'); return;
end
stepSize = round(stepSize);

% check windSize input
if ~isnumeric(windSize)
    warning('Input for ''binSize'' must be numeric. Try again.'); return;
elseif windSize <= 250 || windSize > 5000
    warning('Input for ''binSize'' must be >= 250 or < 5000. Try again.'); return;
end
windSize = round(windSize);

if nargin < 3
    plotOpt = false;
end

%% define year range

yearRange = [100 22000];

if mod(nanmax(yearRange),stepSize) == 0 % trace goes back to 22000, iCESM back to 21000
	ageBin = [nanmin(yearRange):stepSize:nanmax(yearRange)]';
else
	oldYear = nanmax(yearRange) + (binSize - mod(nanmax(yearRange),stepSize));
	ageBin = [nanmin(yearRange):stepSize:stepSize]';
end

% define start and end points within ageBin
[~, indStart] = nanmin(abs(windSize/2 - ageBin)); indStart = indStart+1;
[~, indEnd] = nanmin(abs(nanmax(yearRange - windSize/2) - ageBin)); indEnd = indEnd-1;

%% load proxies;
cd sst_proxy_repository/
    load d18o_proxies_sorted_21-Jan-2020.mat % loads struct 'd18o'
    load mgca_proxies_sorted_21-Jan-2020.mat % loads struct 'mgca'
    load tex86_proxies_sorted_21-Jan-2020.mat% loads struct 'tex86'
    load uk37_proxies_sorted_21-Jan-2020.mat % loads struct 'uk37'
cd ../

% go through each dataset, and remove values outside the yearRange
for i = 1:length(d18o.age)
    indexer = d18o.age{i} < nanmin(yearRange) | d18o.age{i} > nanmax(yearRange);
    d18o.age{i}(indexer) = [];
end
for i = 1:length(mgca.age)
    indexer = mgca.age{i} < nanmin(yearRange) | mgca.age{i} > nanmax(yearRange);
    mgca.age{i}(indexer) = [];
end
for i = 1:length(uk37.age)
    indexer = uk37.age{i} < nanmin(yearRange) | uk37.age{i} > nanmax(yearRange);
    uk37.age{i}(indexer) = [];
end
for i = 1:length(tex86.age)
    indexer = tex86.age{i} < nanmin(yearRange) | tex86.age{i} > nanmax(yearRange);
    tex86.age{i}(indexer) = [];
end

%% Now, run loop
% nbins = 100;
nbins = 0:20:3000;

% d18O
d18o_diffs = [];
d18o_dist = cell(size(ageBin));
for j = indStart:indEnd
    for i = 1:length(d18o.age)
        indexer = d18o.age{i} >= (ageBin(j) - windSize/2) & d18o.age{i} < (ageBin(j) + windSize/2);
        curr = d18o.age{i}(indexer);
        if length(curr) >= 3
            plusdiff = abs(curr(2:end-1) - curr(1:end-2));
            minusdiff = abs(curr(3:end) - curr(2:end-1));
            meandiff = nanmean([plusdiff,minusdiff],2);
            d18o_diffs = [d18o_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
        end
    end
    [d18o_dist{j}.N, d18o_dist{j}.edges] = histcounts(d18o_diffs,nbins,'normalization','cdf'); % cdf pdf
    d18o_dist{j}.edges = (d18o_dist{j}.edges(2:end) + d18o_dist{j}.edges(1:end-1))./2;
    d18o_dist{j}.dist = d18o_diffs;
    d18o_dist{j}.perc90 = prctile(d18o_diffs,90);
    d18o_dist{j}.perc50 = prctile(d18o_diffs,50);
	d18o_diffs = []; % reset
end

% MgCa
mgca_diffs = [];
mgca_dist = cell(size(ageBin));
for j = indStart:indEnd
    for i = 1:length(mgca.age)
        indexer = mgca.age{i} >= (ageBin(j) - windSize/2) & mgca.age{i} < (ageBin(j) + windSize/2);
        curr = mgca.age{i}(indexer);
        if length(curr) >= 3
            plusdiff = abs(curr(2:end-1) - curr(1:end-2));
            minusdiff = abs(curr(3:end) - curr(2:end-1));
            meandiff = nanmean([plusdiff,minusdiff],2);
            mgca_diffs = [mgca_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
        end
    end
    [mgca_dist{j}.N, mgca_dist{j}.edges] = histcounts(mgca_diffs,nbins,'normalization','cdf'); % cdf pdf
    mgca_dist{j}.edges = (mgca_dist{j}.edges(2:end) + mgca_dist{j}.edges(1:end-1))./2;
    mgca_dist{j}.dist = mgca_diffs;
    mgca_dist{j}.perc90 = prctile(mgca_diffs,90);
    mgca_dist{j}.perc50 = prctile(mgca_diffs,50);
	mgca_diffs = []; % reset
end

% Uk37
uk37_diffs = [];
uk37_dist = cell(size(ageBin));
for j = indStart:indEnd
    for i = 1:length(uk37.age)
        indexer = uk37.age{i} >= (ageBin(j) - windSize/2) & uk37.age{i} < (ageBin(j) + windSize/2);
        curr = uk37.age{i}(indexer);
        if length(curr) >= 3
            plusdiff = abs(curr(2:end-1) - curr(1:end-2));
            minusdiff = abs(curr(3:end) - curr(2:end-1));
            meandiff = nanmean([plusdiff,minusdiff],2);
            uk37_diffs = [uk37_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
        end
    end
    [uk37_dist{j}.N, uk37_dist{j}.edges] = histcounts(uk37_diffs,nbins,'normalization','cdf'); % cdf pdf
    uk37_dist{j}.edges = (uk37_dist{j}.edges(2:end) + uk37_dist{j}.edges(1:end-1))./2;
    uk37_dist{j}.dist = uk37_diffs;
    uk37_dist{j}.perc90 = prctile(uk37_diffs,90);
    uk37_dist{j}.perc50 = prctile(uk37_diffs,50);
	uk37_diffs = []; % reset
end

% Tex86
tex86_diffs = [];
tex86_dist = cell(size(ageBin));
for j = indStart:indEnd
    for i = 1:length(tex86.age)
        indexer = tex86.age{i} >= (ageBin(j) - windSize/2) & tex86.age{i} < (ageBin(j) + windSize/2);
        curr = tex86.age{i}(indexer);
        if length(curr) >= 3
            plusdiff = abs(curr(2:end-1) - curr(1:end-2));
            minusdiff = abs(curr(3:end) - curr(2:end-1));
            meandiff = nanmean([plusdiff,minusdiff],2);
            tex86_diffs = [tex86_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
        end
    end
    [tex86_dist{j}.N, tex86_dist{j}.edges] = histcounts(tex86_diffs,nbins,'normalization','cdf'); % cdf pdf
    tex86_dist{j}.edges = (tex86_dist{j}.edges(2:end) + tex86_dist{j}.edges(1:end-1))./2;
    tex86_dist{j}.dist = tex86_diffs;
    tex86_dist{j}.perc90 = prctile(tex86_diffs,90);
    tex86_dist{j}.perc50 = prctile(tex86_diffs,50);
	tex86_diffs = []; % reset
end
   
%% plotting option

cumSum_d18o = nan(length(ageBin),length(nbins)-1); perc_d18o = nan(length(ageBin),2);
cumSum_mgca = nan(length(ageBin),length(nbins)-1); perc_mgca = nan(length(ageBin),2);
cumSum_uk37 = nan(length(ageBin),length(nbins)-1); perc_uk37 = nan(length(ageBin),2);
cumSum_tex86 = nan(length(ageBin),length(nbins)-1);perc_tex86 = nan(length(ageBin),2);
ageBin_plot = nan(size(ageBin)); % ageBin_plot(ageBin < ageBin(indStart) | ageBin > ageBin(indEnd)) = nan;

for j = indStart:indEnd
    cumSum_d18o(j,:) = d18o_dist{j}.N; perc_d18o(j,1:2) = [d18o_dist{j}.perc50, d18o_dist{j}.perc90];
    cumSum_mgca(j,:) = mgca_dist{j}.N; perc_mgca(j,1:2) = [mgca_dist{j}.perc50, mgca_dist{j}.perc90];
    cumSum_uk37(j,:) = uk37_dist{j}.N; perc_uk37(j,1:2) = [uk37_dist{j}.perc50, uk37_dist{j}.perc90];
    cumSum_tex86(j,:) = tex86_dist{j}.N;  perc_tex86(j,1:2) = [tex86_dist{j}.perc50, tex86_dist{j}.perc90];
    ageBin_plot(j) = ageBin(j);
end
edges = tex86_dist{round(size(cumSum_d18o,1)/2)}.edges;

if plotOpt 
    
    % assign colors
    cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
    warning('off','all')
        CT1 = cbrewer('seq','Reds' ,101); 
        CT2 = cbrewer('seq','Blues' ,101);
        CT3 = cbrewer('seq','Purples' ,101);    
        CT4 = cbrewer('seq','Greens' ,101);    
        cd '/Users/matthewosman/Documents/MATLAB/colornames';
        [~,gold] = colornames('CSS','Goldenrod');
    warning('on','all')
    cd(home_fold)

    h = figure; hold on;
        set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
        h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.45*Pix_SS(3),.80*Pix_SS(4)];
        ax = gca; ax.Visible = 'off';

    % d18O
    ax1 = axes('Position',[0.10 0.55 0.28 0.35]); hold on;
        set(ax1,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
        p = pcolor(edges,ageBin-stepSize/2,cumSum_d18o); shading flat; 
        colormap(ax1,CT1); 
        c = colorbar('eastoutside'); c.Label.String = [{'CDF (~)'}]; set(c,'Fontsize',11);
        title('\delta^{18}O','Fontsize',11,'Fontweight','normal'); caxis([0 1]);
        set(gca,'Xscale','log'); xlim([10 1000]); 
            set(gca,'Xtick',[10 50 100 200 400 600 1000]);
        l1 = stairs(perc_d18o(:,1),ageBin,'color',gold,'LineStyle','--','Linewidth',1.5);
        l2 = stairs(perc_d18o(:,2),ageBin,'color',gold,'LineStyle','-','Linewidth',1.5);
        xlabel('Sample Resolution (years)'); ylabel('Years BP');
        set(gca,'Xscale','log'); xlim([10 500]); ylim([0 nanmax(yearRange)])
            set(gca,'Xtick',[10 50 100 200 500]); % 'Ytick',[2500:2500:20000],'YticklabelRotation',45);
        legend([l1 l2],'50^{th}%','90^{th}%','Box','off','Orientation','vertical','Location','southeast');
    
    % mgca
    ax2 = axes('Position',[0.575 0.55 0.28 0.35]); hold on;
        set(ax2,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
        p = pcolor(edges,ageBin-stepSize/2,cumSum_mgca); shading flat; 
        colormap(ax2,CT2); 
        c = colorbar('eastoutside'); c.Label.String = [{'CDF (~)'}]; set(c,'Fontsize',11);
        title('Mg/Ca','Fontsize',11,'Fontweight','normal'); caxis([0 1]);
        set(gca,'Xscale','log'); xlim([10 1000]); 
            set(gca,'Xtick',[10 50 100 200 400 600 1000]);
        l1 = stairs(perc_mgca(:,1),ageBin,'color',gold,'LineStyle','--','Linewidth',1.5);
        l2 = stairs(perc_mgca(:,2),ageBin,'color',gold,'LineStyle','-','Linewidth',1.5);
        xlabel('Sample Resolution (years)'); ylabel('Years BP');
        set(gca,'Xscale','log'); xlim([10 500]); ylim([0 nanmax(yearRange)])
            set(gca,'Xtick',[10 50 100 200 500]); % 'Ytick',[2500:2500:20000],'YticklabelRotation',45);
        legend([l1 l2],'50^{th}%','90^{th}%','Box','off','Orientation','vertical','Location','southeast');
        
    % uk37   
    ax3 = axes('Position',[0.10 0.10 0.28 0.35]); hold on;
        set(ax3,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
        p = pcolor(edges,ageBin-stepSize/2,cumSum_uk37); shading flat; 
        colormap(ax3,CT3); 
        c = colorbar('eastoutside'); c.Label.String = [{'CDF (~)'}]; set(c,'Fontsize',11);
        title('U^{K''}_{37}','Fontsize',11,'Fontweight','normal'); caxis([0 1]);
        set(gca,'Xscale','log'); xlim([10 1000]); 
            set(gca,'Xtick',[10 50 100 200 400 600 1000]);
        l1 = stairs(perc_uk37(:,1),ageBin,'color',gold,'LineStyle','--','Linewidth',1.5);
        l2 = stairs(perc_uk37(:,2),ageBin,'color',gold,'LineStyle','-','Linewidth',1.5);
        xlabel('Sample Resolution (years)'); ylabel('Years BP');
        set(gca,'Xscale','log'); xlim([10 500]); ylim([0 nanmax(yearRange)])
            set(gca,'Xtick',[10 50 100 200 500]); % 'Ytick',[2500:2500:20000],'YticklabelRotation',45);
        legend([l1 l2],'50^{th}%','90^{th}%','Box','off','Orientation','vertical','Location','southeast');
            
    % tex86
    ax4 = axes('Position',[0.575 0.10 0.28 0.35]); hold on;
        set(ax4,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
        p = pcolor(edges,ageBin-stepSize/2,cumSum_tex86); shading flat; 
        colormap(ax4,CT4); 
        c = colorbar('eastoutside'); c.Label.String = [{'CDF (~)'}]; set(c,'Fontsize',11);
        title('TEX^{86}','Fontsize',11,'Fontweight','normal'); caxis([0 1]);
        set(gca,'Xscale','log'); xlim([10 1000]); 
            set(gca,'Xtick',[10 50 100 200 400 600 1000]);
        l1 = stairs(perc_tex86(:,1),ageBin,'color',gold,'LineStyle','--','Linewidth',1.5);
        l2 = stairs(perc_tex86(:,2),ageBin,'color',gold,'LineStyle','-','Linewidth',1.5);
        xlabel('Sample Resolution (years)'); ylabel('Years BP');
        set(gca,'Xscale','log'); xlim([10 500]); ylim([0 nanmax(yearRange)])
            set(gca,'Xtick',[10 50 100 200 500]); % 'Ytick',[2500:2500:20000],'YticklabelRotation',45);
        legend([l1 l2],'50^{th}%','90^{th}%','Box','off','Orientation','vertical','Location','southeast');
        
end
