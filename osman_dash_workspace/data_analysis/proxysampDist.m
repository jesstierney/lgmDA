% Written by M. Osman, Jan2020 (mattosman@arizona.edu) in Matlab2018b
% This function determines the distribution of sampling frequencies (i.e.,
% resolutions) of the LGM - present DA proxy compilation (see Tierney et al., 
% in review).  Insight into the resolution (including potential resolution
% time-dependency) is an important consideration for determining the
% binning size of proxies prior to assimilating
% 
% This function requires:
% 1. A folder named ''sst_proxy_repository/'', housing the following files:
%   a. d18o_proxies_sorted_21-Jan-2020.mat
%   b. mgca_proxies_sorted_21-Jan-2020.mat
%   c. tex86_proxies_sorted_21-Jan-2020.mat
%   d. uk37_proxies_sorted_21-Jan-2020.mat
% NOTE these files can be created via extractProx.m"
clear all; 
home_fold = cd;

%% load proxies;
cd sst_proxy_repository/
    load d18o_proxies_sorted_21-Jan-2020.mat % loads struct 'd18o'
    load mgca_proxies_sorted_21-Jan-2020.mat % loads struct 'mgca'
    load tex86_proxies_sorted_21-Jan-2020.mat% loads struct 'tex86'
    load uk37_proxies_sorted_21-Jan-2020.mat % loads struct 'uk37'
cd ../

% go through each dataset, and remove values outside the yearRange
yearRange = [100 22000]; % years BP, so input of 100 removes all non-preIndust data
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
    
% now, initialize script - 
% important to note that the data are not equally-spaced; to circumvent, 
% we'll take the a = n-1 difference, b = n+1 difference, then mean(a+b) 
d18o_diffs = [];
for i = 1:length(d18o.age)
    plusdiff = abs(d18o.age{i}(2:end-1) - d18o.age{i}(1:end-2));
    minusdiff = abs(d18o.age{i}(3:end) - d18o.age{i}(2:end-1));
    meandiff = nanmean([plusdiff,minusdiff],2);
    d18o_diffs = [d18o_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
end

mgca_diffs = [];
for i = 1:length(mgca.age)
    plusdiff = abs(mgca.age{i}(2:end-1) - mgca.age{i}(1:end-2));
    minusdiff = abs(mgca.age{i}(3:end) - mgca.age{i}(2:end-1));
    meandiff = nanmean([plusdiff,minusdiff],2);
    mgca_diffs = [mgca_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
end

uk37_diffs = [];
for i = 1:length(uk37.age)
    plusdiff = abs(uk37.age{i}(2:end-1) - uk37.age{i}(1:end-2));
    minusdiff = abs(uk37.age{i}(3:end) - uk37.age{i}(2:end-1));
    meandiff = nanmean([plusdiff,minusdiff],2);
    uk37_diffs = [uk37_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
end

tex86_diffs = [];
for i = 1:length(tex86.age)
    plusdiff = abs(tex86.age{i}(2:end-1) - tex86.age{i}(1:end-2));
    minusdiff = abs(tex86.age{i}(3:end) - tex86.age{i}(2:end-1));
    meandiff = nanmean([plusdiff,minusdiff],2);
    tex86_diffs = [tex86_diffs; meandiff]; clearvars plusdiff minusdiff meandiff
end

%% Plot pdf figure;
alphaVal = 0.8;

ncountsmax = 300;
ncount = [round(ncountsmax.*(length(d18o_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(mgca_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(uk37_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(tex86_diffs)/length(d18o_diffs)))];
ncount = [300 300 300 300];
[N.d18o,edges.d18o] = histcounts(d18o_diffs,ncount(1),'normalization','pdf'); % cdf pdf
[N.mgca,edges.mgca] = histcounts(mgca_diffs,ncount(2),'normalization','pdf'); % cdf
[N.uk37,edges.uk37] = histcounts(uk37_diffs,ncount(3),'normalization','pdf'); % cdf
[N.tex86,edges.tex86] = histcounts(tex86_diffs,ncount(4),'normalization','pdf'); % cdf

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,5); 
    CT2 = cbrewer('seq','Blues' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);    
    cd '/Users/matthewosman/Documents/MATLAB/colornames';
    [~,gold] = colornames('CSS','Goldenrod');
warning('on','all')
cd(home_fold)

h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.40*Pix_SS(3),.40*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';

ax1 = axes('Position',[0.15 0.15 0.70 0.70]); hold on;
	set(ax1,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
    b1 = bar((edges.d18o(2:end) + edges.d18o(1:end-1))./2,N.d18o,'Edgecolor','none','FaceColor',CT1(end-1,:),'BarWidth',1,'FaceAlpha',alphaVal); % 0.7.*[1 1 1]
    b2 = bar((edges.mgca(2:end) + edges.mgca(1:end-1))./2,N.mgca,'Edgecolor','none','FaceColor',CT2(end-1,:),'BarWidth',1,'FaceAlpha',alphaVal);
    b3 = bar((edges.uk37(2:end) + edges.uk37(1:end-1))./2,N.uk37,'Edgecolor','none','FaceColor',CT3(end-1,:),'BarWidth',1,'FaceAlpha',alphaVal);
    b4 = bar((edges.tex86(2:end) + edges.tex86(1:end-1))./2,N.tex86,'Edgecolor','none','FaceColor',gold,'BarWidth',1,'FaceAlpha',alphaVal);
%     b1 = plot((edges.d18o(2:end) + edges.d18o(1:end-1))./2,N.d18o,'Linewidth',2,'Color',CT1(end-1,:));
%     b1 = plot((edges.mgca(2:end) + edges.mgca(1:end-1))./2,N.mgca,'Linewidth',2,'Color',CT2(end-1,:));
%     b1 = plot((edges.uk37(2:end) + edges.uk37(1:end-1))./2,N.uk37,'Linewidth',2,'Color',CT3(end-1,:));
%     b1 = plot((edges.tex86(2:end) + edges.tex86(1:end-1))./2,N.tex86,'Linewidth',2,'Color',gold);
    xlabel('Sample Resolution (years)'); ylabel('PDF (~)');
    set(gca,'Xscale','log'); xlim([10 1000]); 
    set(gca,'Xtick',[10 50 100 200 400 600 1000]);
    % set(gca,'Xscale','linear'); xlim([0 1000]); 
    legend([b1 b2 b3 b4],'\delta^{18}O','Mg/Ca','U^{K''}_{37}','TEX^{86}','Box','off','Orientation','vertical','Location','northeast');
    
%% Plot cdf figure;
alphaVal = 0.9;

ncountsmax = 300;
ncount = [round(ncountsmax.*(length(d18o_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(mgca_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(uk37_diffs)/length(d18o_diffs))),...
          round(ncountsmax.*(length(tex86_diffs)/length(d18o_diffs)))];
ncount = [300 300 300 300];
[N.d18o,edges.d18o] = histcounts(d18o_diffs,ncount(1),'normalization','cdf'); % cdf pdf
[N.mgca,edges.mgca] = histcounts(mgca_diffs,ncount(2),'normalization','cdf'); % cdf
[N.uk37,edges.uk37] = histcounts(uk37_diffs,ncount(3),'normalization','cdf'); % cdf
[N.tex86,edges.tex86] = histcounts(tex86_diffs,ncount(4),'normalization','cdf'); % cdf

% assign colors
cd '/Users/matthewosman/Documents/MATLAB/cbrewer' 
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,5); 
    CT2 = cbrewer('seq','Blues' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);    
    cd '/Users/matthewosman/Documents/MATLAB/colornames';
    [~,gold] = colornames('CSS','Goldenrod');
warning('on','all')
cd(home_fold)

h = figure; hold on;
    set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
    h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.40*Pix_SS(3),.40*Pix_SS(4)];
    ax = gca; ax.Visible = 'off';

ax1 = axes('Position',[0.15 0.15 0.70 0.70]); hold on;
	set(ax1,'Yaxislocation','left','xaxisLocation','bottom','fontsize',11,'box','on','Color','None','TickDir','in','linewidth',1.5); %'Xcolor','none','Xtick',[])
    l1 = line([0.1 1000],[0.9 0.9],'color',0.3.*[1 1 1],'LineStyle','--','Linewidth',1.5);
    p1 = plot((edges.d18o(2:end) + edges.d18o(1:end-1))./2,N.d18o,'Linewidth',2,'Color',CT1(end-1,:));
    p2 = plot((edges.mgca(2:end) + edges.mgca(1:end-1))./2,N.mgca,'Linewidth',2,'Color',CT2(end-1,:));
    p3 = plot((edges.uk37(2:end) + edges.uk37(1:end-1))./2,N.uk37,'Linewidth',2,'Color',CT3(end-1,:));
    p4 = plot((edges.tex86(2:end) + edges.tex86(1:end-1))./2,N.tex86,'Linewidth',2,'Color',gold);
    xlabel('Sample Resolution (years)'); ylabel('CDF (~)');
    set(gca,'Xscale','log'); xlim([10 1000]); 
        set(gca,'Xtick',[10 50 100 200 400 600 1000]);
    % set(gca,'Xscale','linear'); xlim([0 1000]); 
    legend([p1 p2 p3 p4],'\delta^{18}O','Mg/Ca','U^{K''}_{37}','TEX^{86}','Box','off','Orientation','vertical','Location','southeast');