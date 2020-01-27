% This script extracts all the proxy data from the .nc file, 'proxysite-2019-10-09.nc'
% and stores it into a collection of .mat files; 
% WARNING! this file takes quite a bit of time to run.....
% To run, requires :
%   1. the file "proxysite-2019-10-09.nc"
%   2. a folder named "sst_proxy_repository" to place sorted proxy data into
%   3. access to the BAYMAG package (https://github.com/jesstierney/BAYMAG); 
%       note will need to change lines 304/308
% Written by M. Osman (Nov 2019; mattosman@arizona.edu) via Matlab2018b
clear all; 
home_fold = cd;
%% proxy types: 
% 'uk37', 
uk37_name = {'uk37'};
%
% 'd18o_acicula', 'd18o_bulloides', 'd18o_dutertrei', 'd18o_inflata', 'd18o_mabahethi', 
%   'd18o_marginata', 'd18o_menardii', 'd18o_obliquiloculata', 'd18o_pachyderma', 'd18o_pachyderma_d',
%   'd18o_quinqueloba','d18o_ruber', 'd18o_ruber_pink', 'd18o_sacculifer', 'd18o_tumida',
d18O_name = {'d18o_acicula'; 'd18o_bulloides'; 'd18o_dutertrei'; 'd18o_inflata'; 'd18o_mabahethi'; ...
    'd18o_marginata'; 'd18o_menardii'; 'd18o_obliquiloculata'; 'd18o_pachyderma'; 'd18o_pachyderma_d'; ...
	'd18o_quinqueloba'; 'd18o_ruber'; 'd18o_ruber_pink'; 'd18o_sacculifer'; 'd18o_tumida'};
% 
% 'd13c_bulloides', 'd13c_dutertrei', 'd13c_pachyderma', 'd13c_ruber_pink', 'd13c_sacculifer', 
d13c_name = {'d13c_bulloides'; 'd13c_dutertrei'; 'd13c_pachyderma'; 'd13c_ruber_pink'; 'd13c_sacculifer'};
% 
% 'mgca_bulloides', 'mgca_crassaformis', 'mgca_dutertrei', 'mgca_obliquiloculata', 'mgca_pachyderma_d', 'mgca_ruber', 
%   'mgca_ruber_lato', 'mgca_ruber_stricto', 'mgca_sacculifer', 'mgca_truncatulinoides'
mgca_name = {'mgca_bulloides'; 'mgca_crassaformis'; 'mgca_dutertrei'; 'mgca_obliquiloculata'; 'mgca_pachyderma_d';... 
    'mgca_ruber'; 'mgca_ruber_lato'; 'mgca_ruber_stricto'; 'mgca_sacculifer'; 'mgca_truncatulinoides'};
%  
% 'tex86'
tex86_name = {'tex86'};
% 
% 'percent_pachysin'
percent_pachysin_name = {'percent_pachysin'};

%% Extract each of the variable types
filename='proxysite-2019-10-09.nc';
finfo = ncinfo(filename);
cores = {finfo.Groups.Name}; % get core names
        
% mgca
disp('Starting ''mgca''')
h = waitbar(0,'Please wait... sites remaining:');
identifier = 'mgca';
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ mgca
        if nansum(indexer) > 0 % if mgca exists within the data
            indices = find(indexer); % find indices of the mgca variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                mgca.coreName{m,1}      = cores{i};
                mgca.varName{m,1}       = curr_name;
                mgca.age{m,1}           = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                mgca.foram_species{m,1} = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'foraminifera_type');
                mgca.cleaning{m,1}      = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'mgca_cleaning_protocol');
                mgca.units{m,1}         = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'units');                                
                mgca.data{m,1}          = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                mgca.elev(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % mgca.latlon(m,[1,2])  = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                mgca.lats(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                mgca.lons(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end   
close(h)
cd sst_proxy_repository/; save(['mgca_proxies_sorted_',date,'.mat'],'mgca'); cd ../
return; 

% d18o
identifier = 'd18o';
disp('Starting ''d18o''')
h = waitbar(0,'Please wait... sites remaining:');
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ d18o
        if nansum(indexer) > 0 % if d18O exists within the data
            indices = find(indexer); % find indices of the d18O variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                d18o.coreName{m,1}      = cores{i};
                d18o.varName{m,1}       = curr_name;
                d18o.age{m,1}           = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                d18o.foram_species{m,1} = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'foraminifera_type');
                d18o.units{m,1}         = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'units');
                d18o.data{m,1}          = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                d18o.elev(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % d18o.latlon(m,[1,2])  = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                d18o.lats(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                d18o.lons(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end
close(h)
cd sst_proxy_repository/; save(['d18o_proxies_sorted_',date,'.mat'],'d18o'); cd ../
    

identifier = 'uk37';
disp('Starting ''uk37''')
h = waitbar(0,'Please wait... sites remaining:');
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ uk37
        if nansum(indexer) > 0 % if uk37 exists within the data
            indices = find(indexer); % find indices of the uk37 variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                uk37.coreName{m,1}   = cores{i};
                uk37.varName{m,1}    = curr_name;
                uk37.age{m,1}        = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                uk37.units{m,1}      = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'units');
                uk37.data{m,1}       = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                uk37.elev(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % uk37.latlon(m,[1,2]) = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                uk37.lats(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                uk37.lons(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end
close(h)
cd sst_proxy_repository/; save(['uk37_proxies_sorted_',date,'.mat'],'uk37'); cd ../
    
      
% tex86
disp('Starting ''tex86''')
identifier = 'tex86';
h = waitbar(0,'Please wait... sites remaining:');
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ tex86
        if nansum(indexer) > 0 % if tex86 exists within the data
            indices = find(indexer); % find indices of the tex86 variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                tex86.coreName{m,1}   = cores{i};
                tex86.varName{m,1}    = curr_name;
                tex86.age{m,1}        = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                tex86.units{m,1}      = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'units');
                tex86.data{m,1}       = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                tex86.elev(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % tex86.latlon(m,[1,2]) = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                tex86.lats(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                tex86.lons(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end 
close(h)
cd sst_proxy_repository/; save(['tex86_proxies_sorted_',date,'.mat'],'tex86'); cd ../
    

% d13c
disp('Starting ''d13c''')
h = waitbar(0,'Please wait... sites remaining:');
identifier = 'd13c';
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ d13c
        if nansum(indexer) > 0 % if d13c exists within the data
            indices = find(indexer); % find indices of the d13c variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                d13c.coreName{m,1}      = cores{i};
                d13c.varName{m,1}       = curr_name;
                d13c.age{m,1}           = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                d13c.foram_species{m,1} = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'foraminifera_type');
                d13c.units{m,1}         = ncreadatt(filename,strcat('/',cores{i},'/data/',curr_name),'units');                
                d13c.data{m,1}          = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                d13c.elev(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % d13c.latlon(m,[1,2])  = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                d13c.lats(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                d13c.lons(m,1)          = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end
close(h)
cd sst_proxy_repository/; save(['d13c_proxies_sorted_',date,'.mat'],'d13c'); cd ../


% percent_pachysin
disp('Starting ''percent_pachysin''')
identifier = 'percent_pachysin';
h = waitbar(0,'Please wait... sites remaining:');
    m = 1; % inititiate
    for i = 1:length(cores) % for each site
        waitbar(i / length(cores));
        indexer = startsWith({finfo.Groups(i).Groups(2).Variables.Name},identifier); % identifies which, if any, data strings start w/ percent_pachysin
        if nansum(indexer) > 0 % if percent_pachysin exists within the data
            indices = find(indexer); % find indices of the percent_pachysin variables
            for j = 1:length(indices)
                curr_names = {finfo.Groups(i).Groups(2).Variables.Name};
                curr_name  = curr_names{indices(j)};
                percent_pachysin.coreName{m,1}   = cores{i};
                percent_pachysin.varName{m,1}    = curr_name;
                percent_pachysin.age{m,1}        = ncread(filename,strcat('/',cores{i},'/data/age_median'));
                percent_pachysin.data{m,1}       = ncread(filename,strcat('/',cores{i},'/data/',curr_name));
                percent_pachysin.elev(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'elevation');
                % percent_pachysin.latlon(m,[1,2]) = [ncreadatt(filename,strcat('/',cores{i},'/'),'latitude'), ncreadatt(filename,strcat('/',cores{i},'/'),'longitude')];
                percent_pachysin.lats(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'latitude');
                percent_pachysin.lons(m,1)       = ncreadatt(filename,strcat('/',cores{i},'/'),'longitude');
                m = m+1; % reset index
                clearvars curr_names curr_name
            end
        end   
    end 
close(h)
cd sst_proxy_repository/; save(['percent_pachysin_proxies_sorted_',date,'.mat'],'percent_pachysin'); cd ../


%% Additional variables needed for PSM's

cd('sst_proxy_repository/');
    load mgca_proxies_sorted_21-Jan-2020.mat
    load d18o_proxies_sorted_21-Jan-2020.mat
cd ../

% Species specifications for d18o and MgCa:

% from baymag / bayspar; input parameter:
% species = string of target species. six options:
%   'ruber' = Globigerinoides ruber, white or pink
%   'bulloides' = Globigerina bulloides
%   'sacculifer' = Trilobatus sacculifer
%   'pachy' = Neogloboquadrina pachyderma
%   'incompta' = Neogloboquadrina incompta
%   'all' = pooled calibration, annual SST
%   'all_sea' = pooled calibration, seasonal SST

    for i = 1:length(mgca.varName)
        % only need to extract string elements 6-10 from "varName" to identify species
        currString = mgca.varName{i}(6:10);
        if strcmp(currString,'ruber')
            mgca.species{i,1} = 'ruber';
        elseif strcmp(currString,'bullo')
            mgca.species{i,1} = 'bulloides';
        elseif strcmp(currString,'saccu')
            mgca.species{i,1} = 'sacculifer';
        elseif strcmp(currString,'pachy')
            mgca.species{i,1} = 'pachy';
        elseif strcmp(currString,'incom')
            mgca.species{i,1} = 'incompta';
        else
            mgca.species{i,1} = 'other';
        end
    end

    for i = 1:length(d18o.varName)
        % only need to extract string elements 6-10 from "varName" to identify species
        currString = d18o.varName{i}(6:10);
        if strcmp(currString,'ruber')
            d18o.species{i,1} = 'ruber';
        elseif strcmp(currString,'bullo')
            d18o.species{i,1} = 'bulloides';
        elseif strcmp(currString,'saccu')
            d18o.species{i,1} = 'sacculifer';
        elseif strcmp(currString,'pachy')
            d18o.species{i,1} = 'pachy';
        elseif strcmp(currString,'incom')
            d18o.species{i,1} = 'incompta';
        else
            d18o.species{i,1} = 'other';
        end
    end

% Now, for mgca, need to determine cleaning method; from baymag:
% clean = scalar to describe cleaning technique. options:
%   1 = reductive 
%   0 = oxidative
    for i = 1:length(mgca.varName)
        currString = mgca.cleaning{i};
        if strcmp(currString,'Barker cleaning with hydrogen peroxide')
            mgca.cleaningTF(i,1) = 0;
        elseif strcmp(currString,'Fully reductive cleaning')
            mgca.cleaningTF(i,1) = 1;
        else
            warning(['No cleanng condition matched for siteID-',num2str(i)]);
        end
    end
    mgca=rmfield(mgca,'cleaning');
    mgca.cleaning = mgca.cleaningTF;
    mgca=rmfield(mgca,'cleaningTF');

% finally, need ph and omega - NOTE - pH still needs to be corrected for sample age !
    cd('/Users/matthewosman/Documents/GitHub/DASH/3. PSMs/Specific Forward Models/BAYMAG-master/tsget')
        for i = 1:length(mgca.lats)
            [mgca.omega(i,1), mgca.ph(i,1)] = omgph(mgca.lats(i),mgca.lons(i),abs(double(mgca.elev(i))));
            % [mgca.omega(i,1), mgca.ph(i,1)] = omgph(mgca.lats(i),mgca.lons(i),abs(double(mgca.elev(i))));
        end
    cd('/Users/matthewosman/Documents/GitHub/DASH');

% resave the mgca and d18o structures!
cd sst_proxy_repository/; save(['mgca_proxies_sorted_',date,'.mat'],'mgca'); cd ../
cd sst_proxy_repository/; save(['d18o_proxies_sorted_',date,'.mat'],'d18o'); cd ../
