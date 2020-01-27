function [allData, calibData, validData] = prepProx(calibFrac,binSize) 
% Written by M. Osman, Jan2020 (mattosman@arizona.edu) via Matlab2018b
% This function takes motivation from J. Tierney's "preProxies.m" function.
%   In brief, this function allocates the proxy data into 3 structures sorted
%   by time period relative to 0 yrs BP, such that the 3 structures contain
%   a) all available data, b) some fraction of this proxy data randomly set aside for 
%   calibraton, and c) the remaining proxy data, set aside for validation:
% INPUTS
% calibFrac : the fraction of data set aside for validation; must be
%   between >= 0 and < 1
% binSize : integer denoting bin size (e.g., input of binSize = 50 would
%   find all data falling within interval -1 to -50 yr BP, -51 to -100 yr
%   BP, -101 to -150 yr BP, etc etc ...); must be >0 and <= 5000

% Requirements to run -- in the Current Folder, you must have:
%   1. a folder named "sst_proxy_repository", containing:
%       a. uk37_proxies_sorted_21-Jan-2020.mat
%       b. tex86_proxies_sorted_21-Jan-2020.mat
%       c. mgca_proxies_sorted_21-Jan-2020.mat
%       d. d18o_proxies_sorted_21-Jan-2020.mat
% NOTE - the above .mat files were in turn created via the function "extract_proxy_data.m" 

% could also add option to search for data in some specific area!
% to do so, would first remove all values from .mat files not within this
% area, before sorting...

% load in the proxy data:
cd sst_proxy_repository
    load uk37_proxies_sorted_21-Jan-2020.mat % uk37
    load tex86_proxies_sorted_21-Jan-2020.mat% tex86
    load d18o_proxies_sorted_21-Jan-2020.mat % d18O
    load mgca_proxies_sorted_21-Jan-2020.mat % mgca
cd ../

if calibFrac < 0 || calibFrac >= 1
    warning('Input for ''calibFrac'' must be >= 0 and < 1! Try again.'); return;
end

if ~isnumeric(binSize)
    warning('Input for ''binSize'' must be numeric. Try again.'); return;
elseif binSize <= 0 || binSize > 5000
    warning('Input for ''binSize'' must be >= 0 or < 5000. Try again.'); return;
end
binSize = round(binSize);

% define year range
if mod(22000,binSize) == 0 % trace goes back to 22000, iCESM back to 21000
	ageBin = [0:binSize:22000]';
else
	oldYear = 22000 + (binSize - mod(22000,binSize));
	ageBin = [0:binSize:oldYear]';
end
    
%% Sort all data:

for i = 1:length(ageBin)-1
    binRange = [ageBin(i)+1, ageBin(i+1)];
    disp(['Sorting proxy samples for ',num2str(nanmin(binRange)),'-',num2str(nanmax(binRange)),' years BP.']);
    % mgca
    allData.mgca.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(mgca.age)
        ind = mgca.age{j} >= ageBin(i)+1 & mgca.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            allData.mgca.age{i,1}(m:(m+indSum-1),1) = mgca.age{j}(ind);
            allData.mgca.data{i,1}(m:(m+indSum-1),1) = mgca.data{j}(ind);
            allData.mgca.elev{i,1}(m:(m+indSum-1),1) = mgca.elev(j);
            allData.mgca.lats{i,1}(m:(m+indSum-1),1) = mgca.lats(j);
            allData.mgca.lons{i,1}(m:(m+indSum-1),1) = mgca.lons(j);
            allData.mgca.cleaning{i,1}(m:(m+indSum-1),1) = mgca.cleaning(j);
            allData.mgca.ph{i,1}(m:(m+indSum-1),1) = mgca.ph(j);
            allData.mgca.omega{i,1}(m:(m+indSum-1),1) = mgca.omega(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.mgca.coreName{i,1}{k,1} = mgca.coreName{j};
                    allData.mgca.species{i,1}{k,1} = mgca.species{j};
                    end
                else
                allData.mgca.coreName{i,1}{m:(m+indSum-1),1} = mgca.coreName{j};
                allData.mgca.species{i,1}{m:(m+indSum-1),1} = mgca.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % d18o
    allData.d18o.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(d18o.age)
        ind = d18o.age{j} >= ageBin(i)+1 & d18o.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            allData.d18o.age{i,1}(m:(m+indSum-1),1) = d18o.age{j}(ind);
            allData.d18o.data{i,1}(m:(m+indSum-1),1) = d18o.data{j}(ind);
            allData.d18o.elev{i,1}(m:(m+indSum-1),1) = d18o.elev(j);
            allData.d18o.lats{i,1}(m:(m+indSum-1),1) = d18o.lats(j);
            allData.d18o.lons{i,1}(m:(m+indSum-1),1) = d18o.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.d18o.coreName{i,1}{k,1} = d18o.coreName{j};
                    allData.d18o.species{i,1}{k,1} = d18o.species{j};
                    end
                else
                allData.d18o.coreName{i,1}{m:(m+indSum-1),1} = d18o.coreName{j};
                allData.d18o.species{i,1}{m:(m+indSum-1),1} = d18o.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % uk37
    allData.uk37.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(uk37.age)
        ind = uk37.age{j} >= ageBin(i)+1 & uk37.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            allData.uk37.age{i,1}(m:(m+indSum-1),1) = uk37.age{j}(ind);
            allData.uk37.data{i,1}(m:(m+indSum-1),1) = uk37.data{j}(ind);
            allData.uk37.elev{i,1}(m:(m+indSum-1),1) = uk37.elev(j);
            allData.uk37.lats{i,1}(m:(m+indSum-1),1) = uk37.lats(j);
            allData.uk37.lons{i,1}(m:(m+indSum-1),1) = uk37.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.uk37.coreName{i,1}{k,1} = uk37.coreName{j};
                    end
                else
                allData.uk37.coreName{i,1}{m:(m+indSum-1),1} = uk37.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end 
    % tex86
    allData.tex86.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(tex86.age)
        ind = tex86.age{j} >= ageBin(i)+1 & tex86.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            allData.tex86.age{i,1}(m:(m+indSum-1),1) = tex86.age{j}(ind);
            allData.tex86.data{i,1}(m:(m+indSum-1),1) = tex86.data{j}(ind);
            allData.tex86.elev{i,1}(m:(m+indSum-1),1) = tex86.elev(j);
            allData.tex86.lats{i,1}(m:(m+indSum-1),1) = tex86.lats(j);
            allData.tex86.lons{i,1}(m:(m+indSum-1),1) = tex86.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.tex86.coreName{i,1}{k,1} = tex86.coreName{j};
                    end
                else
                allData.tex86.coreName{i,1}{m:(m+indSum-1),1} = tex86.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end    
end

%% Now, randomly choose "calibFrac" of original sites for calibration

% mgca
ind_mgca_cal = sort(randperm(length(mgca.age),round(calibFrac.*length(mgca.age))),'ascend');
ind_mgca_val = (1:length(mgca.age)); ind_mgca_val(ind_mgca_cal) = [];
%     ver = unique([(1:length(mgca.age)),ind_mgca_cal]);
% ind_mgca_val = ver(histcounts([(1:length(mgca.age)),ind_mgca_cal],[ver,inf])==1);
    % assign calib data
    mgca_calib.coreName = mgca.coreName(ind_mgca_cal);
    mgca_calib.varName = mgca.varName(ind_mgca_cal);
    mgca_calib.age = mgca.age(ind_mgca_cal);
    mgca_calib.data = mgca.data(ind_mgca_cal);
    mgca_calib.elev = mgca.elev(ind_mgca_cal);
    mgca_calib.lats = mgca.lats(ind_mgca_cal);
    mgca_calib.lons = mgca.lons(ind_mgca_cal);
    mgca_calib.species = mgca.species(ind_mgca_cal);
    mgca_calib.cleaning = mgca.cleaning(ind_mgca_cal);
    mgca_calib.omega = mgca.omega(ind_mgca_cal);
    mgca_calib.ph = mgca.ph(ind_mgca_cal);
    % assign valid data
    mgca_valid.coreName = mgca.coreName(ind_mgca_val);
    mgca_valid.varName = mgca.varName(ind_mgca_val);
    mgca_valid.age = mgca.age(ind_mgca_val);
    mgca_valid.data = mgca.data(ind_mgca_val);
    mgca_valid.elev = mgca.elev(ind_mgca_val);
    mgca_valid.lats = mgca.lats(ind_mgca_val);
    mgca_valid.lons = mgca.lons(ind_mgca_val);
    mgca_valid.species = mgca.species(ind_mgca_val);
    mgca_valid.cleaning = mgca.cleaning(ind_mgca_val);
    mgca_valid.omega = mgca.omega(ind_mgca_val);
    mgca_valid.ph = mgca.ph(ind_mgca_val);

% d18o
ind_d18o_cal = sort(randperm(length(d18o.age),round(calibFrac.*length(d18o.age))),'ascend');
ind_d18o_val = (1:length(d18o.age)); ind_d18o_val(ind_d18o_cal) = [];
%     ver = unique([(1:length(d18o.age)),ind_d18o_cal]);
% ind_d18o_val = ver(histcounts([(1:length(d18o.age)),ind_d18o_cal],[ver,inf])==1);
    % assign calib data
    d18o_calib.coreName = d18o.coreName(ind_d18o_cal);
    d18o_calib.varName = d18o.varName(ind_d18o_cal);
    d18o_calib.age = d18o.age(ind_d18o_cal);
    d18o_calib.data = d18o.data(ind_d18o_cal);
    d18o_calib.elev = d18o.elev(ind_d18o_cal);
    d18o_calib.lats = d18o.lats(ind_d18o_cal);
    d18o_calib.lons = d18o.lons(ind_d18o_cal);
    d18o_calib.species = d18o.species(ind_d18o_cal);
    % assign valid data
    d18o_valid.coreName = d18o.coreName(ind_d18o_val);
    d18o_valid.varName = d18o.varName(ind_d18o_val);
    d18o_valid.age = d18o.age(ind_d18o_val);
    d18o_valid.data = d18o.data(ind_d18o_val);
    d18o_valid.elev = d18o.elev(ind_d18o_val);
    d18o_valid.lats = d18o.lats(ind_d18o_val);
    d18o_valid.lons = d18o.lons(ind_d18o_val);
    d18o_valid.species = d18o.species(ind_d18o_val);

% tex86
ind_tex86_cal = sort(randperm(length(tex86.age),round(calibFrac.*length(tex86.age))),'ascend');
ind_tex86_val = (1:length(tex86.age)); ind_tex86_val(ind_tex86_cal) = [];
%     ver = unique([(1:length(tex86.age)),ind_tex86_cal]);
% ind_tex86_val = ver(histcounts([(1:length(tex86.age)),ind_tex86_cal],[ver,inf])==1);
    % assign calib data
    tex86_calib.coreName = tex86.coreName(ind_tex86_cal);
    tex86_calib.varName = tex86.varName(ind_tex86_cal);
    tex86_calib.age = tex86.age(ind_tex86_cal);
    tex86_calib.data = tex86.data(ind_tex86_cal);
    tex86_calib.elev = tex86.elev(ind_tex86_cal);
    tex86_calib.lats = tex86.lats(ind_tex86_cal);
    tex86_calib.lons = tex86.lons(ind_tex86_cal);
    % assign valid data
    tex86_valid.coreName = tex86.coreName(ind_tex86_val);
    tex86_valid.varName = tex86.varName(ind_tex86_val);
    tex86_valid.age = tex86.age(ind_tex86_val);
    tex86_valid.data = tex86.data(ind_tex86_val);
    tex86_valid.elev = tex86.elev(ind_tex86_val);
    tex86_valid.lats = tex86.lats(ind_tex86_val);
    tex86_valid.lons = tex86.lons(ind_tex86_val);

% uk37
ind_uk37_cal = sort(randperm(length(uk37.age),round(calibFrac.*length(uk37.age))),'ascend');
ind_uk37_val = (1:length(uk37.age)); ind_uk37_val(ind_uk37_cal) = [];
%     ver = unique([(1:length(uk37.age)),ind_uk37_cal]);
% ind_uk37_val = ver(histcounts([(1:length(uk37.age)),ind_uk37_cal],[ver,inf])==1);
    % assign calib data
    uk37_calib.coreName = uk37.coreName(ind_uk37_cal);
    uk37_calib.varName = uk37.varName(ind_uk37_cal);
    uk37_calib.age = uk37.age(ind_uk37_cal);
    uk37_calib.data = uk37.data(ind_uk37_cal);
    uk37_calib.elev = uk37.elev(ind_uk37_cal);
    uk37_calib.lats = uk37.lats(ind_uk37_cal);
    uk37_calib.lons = uk37.lons(ind_uk37_cal);
    % assign valid data
    uk37_valid.coreName = uk37.coreName(ind_uk37_val);
    uk37_valid.varName = uk37.varName(ind_uk37_val);
    uk37_valid.age = uk37.age(ind_uk37_val);
    uk37_valid.data = uk37.data(ind_uk37_val);
    uk37_valid.elev = uk37.elev(ind_uk37_val);
    uk37_valid.lats = uk37.lats(ind_uk37_val);
    uk37_valid.lons = uk37.lons(ind_uk37_val);

%% Sort calib data

for i = 1:length(ageBin)-1
    binRange = [ageBin(i)+1, ageBin(i+1)];
    disp(['Sorting calibration proxy samples for ',num2str(nanmin(binRange)),'-',num2str(nanmax(binRange)),' years BP.']);
    % mgca
    calibData.mgca.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(mgca_calib.age)
        ind = mgca_calib.age{j} >= ageBin(i)+1 & mgca_calib.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            calibData.mgca.age{i,1}(m:(m+indSum-1),1) = mgca_calib.age{j}(ind);
            calibData.mgca.data{i,1}(m:(m+indSum-1),1) = mgca_calib.data{j}(ind);
            calibData.mgca.elev{i,1}(m:(m+indSum-1),1) = mgca_calib.elev(j);
            calibData.mgca.lats{i,1}(m:(m+indSum-1),1) = mgca_calib.lats(j);
            calibData.mgca.lons{i,1}(m:(m+indSum-1),1) = mgca_calib.lons(j);
            calibData.mgca.cleaning{i,1}(m:(m+indSum-1),1) = mgca_calib.cleaning(j);
            calibData.mgca.ph{i,1}(m:(m+indSum-1),1) = mgca_calib.ph(j);
            calibData.mgca.omega{i,1}(m:(m+indSum-1),1) = mgca_calib.omega(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    calibData.mgca.coreName{i,1}{k,1} = mgca_calib.coreName{j};
                    calibData.mgca.species{i,1}{k,1} = mgca_calib.species{j};
                    end
                else
                calibData.mgca.coreName{i,1}{m:(m+indSum-1),1} = mgca_calib.coreName{j};
                calibData.mgca.species{i,1}{m:(m+indSum-1),1} = mgca_calib.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % d18o
    calibData.d18o.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(d18o_calib.age)
        ind = d18o_calib.age{j} >= ageBin(i)+1 & d18o_calib.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            calibData.d18o.age{i,1}(m:(m+indSum-1),1) = d18o_calib.age{j}(ind);
            calibData.d18o.data{i,1}(m:(m+indSum-1),1) = d18o_calib.data{j}(ind);
            calibData.d18o.elev{i,1}(m:(m+indSum-1),1) = d18o_calib.elev(j);
            calibData.d18o.lats{i,1}(m:(m+indSum-1),1) = d18o_calib.lats(j);
            calibData.d18o.lons{i,1}(m:(m+indSum-1),1) = d18o_calib.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.d18o.coreName{i,1}{k,1} = d18o_calib.coreName{j};
                    allData.d18o.species{i,1}{k,1} = d18o_calib.species{j};
                    end
                else
                allData.d18o.coreName{i,1}{m:(m+indSum-1),1} = d18o_calib.coreName{j};
                allData.d18o.species{i,1}{m:(m+indSum-1),1} = d18o_calib.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % uk37
    allData.uk37.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(uk37_calib.age)
        ind = uk37_calib.age{j} >= ageBin(i)+1 & uk37_calib.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            calibData.uk37.age{i,1}(m:(m+indSum-1),1) = uk37_calib.age{j}(ind);
            calibData.uk37.data{i,1}(m:(m+indSum-1),1) = uk37_calib.data{j}(ind);
            calibData.uk37.elev{i,1}(m:(m+indSum-1),1) = uk37_calib.elev(j);
            calibData.uk37.lats{i,1}(m:(m+indSum-1),1) = uk37_calib.lats(j);
            calibData.uk37.lons{i,1}(m:(m+indSum-1),1) = uk37_calib.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.uk37.coreName{i,1}{k,1} = uk37_calib.coreName{j};
                    end
                else
                calibData.uk37.coreName{i,1}{m:(m+indSum-1),1} = uk37_calib.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end 
    % tex86
    calibData.tex86.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(tex86_calib.age)
        ind = tex86_calib.age{j} >= ageBin(i)+1 & tex86_calib.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            calibData.tex86.age{i,1}(m:(m+indSum-1),1) = tex86_calib.age{j}(ind);
            calibData.tex86.data{i,1}(m:(m+indSum-1),1) = tex86_calib.data{j}(ind);
            calibData.tex86.elev{i,1}(m:(m+indSum-1),1) = tex86_calib.elev(j);
            calibData.tex86.lats{i,1}(m:(m+indSum-1),1) = tex86_calib.lats(j);
            calibData.tex86.lons{i,1}(m:(m+indSum-1),1) = tex86_calib.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    calibData.tex86.coreName{i,1}{k,1} = tex86_calib.coreName{j};
                    end
                else
                calibData.tex86.coreName{i,1}{m:(m+indSum-1),1} = tex86_calib.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end    
end

%% Sort valid data

for i = 1:length(ageBin)-1
    binRange = [ageBin(i)+1, ageBin(i+1)];
    disp(['Sorting validation proxy samples for ',num2str(nanmin(binRange)),'-',num2str(nanmax(binRange)),' years BP.']);
    % mgca
    validData.mgca.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(mgca_valid.age)
        ind = mgca_valid.age{j} >= ageBin(i)+1 & mgca_valid.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            validData.mgca.age{i,1}(m:(m+indSum-1),1) = mgca_valid.age{j}(ind);
            validData.mgca.data{i,1}(m:(m+indSum-1),1) = mgca_valid.data{j}(ind);
            validData.mgca.elev{i,1}(m:(m+indSum-1),1) = mgca_valid.elev(j);
            validData.mgca.lats{i,1}(m:(m+indSum-1),1) = mgca_valid.lats(j);
            validData.mgca.lons{i,1}(m:(m+indSum-1),1) = mgca_valid.lons(j);
            validData.mgca.cleaning{i,1}(m:(m+indSum-1),1) = mgca_valid.cleaning(j);
            validData.mgca.ph{i,1}(m:(m+indSum-1),1) = mgca_valid.ph(j);
            validData.mgca.omega{i,1}(m:(m+indSum-1),1) = mgca_valid.omega(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    validData.mgca.coreName{i,1}{k,1} = mgca_valid.coreName{j};
                    validData.mgca.species{i,1}{k,1} = mgca_valid.species{j};
                    end
                else
                validData.mgca.coreName{i,1}{m:(m+indSum-1),1} = mgca_valid.coreName{j};
                validData.mgca.species{i,1}{m:(m+indSum-1),1} = mgca_valid.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % d18o
    validData.d18o.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(d18o_valid.age)
        ind = d18o_valid.age{j} >= ageBin(i)+1 & d18o_valid.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            validData.d18o.age{i,1}(m:(m+indSum-1),1) = d18o_valid.age{j}(ind);
            validData.d18o.data{i,1}(m:(m+indSum-1),1) = d18o_valid.data{j}(ind);
            validData.d18o.elev{i,1}(m:(m+indSum-1),1) = d18o_valid.elev(j);
            validData.d18o.lats{i,1}(m:(m+indSum-1),1) = d18o_valid.lats(j);
            validData.d18o.lons{i,1}(m:(m+indSum-1),1) = d18o_valid.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.d18o.coreName{i,1}{k,1} = d18o_valid.coreName{j};
                    allData.d18o.species{i,1}{k,1} = d18o_valid.species{j};
                    end
                else
                allData.d18o.coreName{i,1}{m:(m+indSum-1),1} = d18o_valid.coreName{j};
                allData.d18o.species{i,1}{m:(m+indSum-1),1} = d18o_valid.species{j};
                end
            % update m
            m = m+indSum-1;
        end
    end
    % uk37
    allData.uk37.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(uk37_valid.age)
        ind = uk37_valid.age{j} >= ageBin(i)+1 & uk37_valid.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            validData.uk37.age{i,1}(m:(m+indSum-1),1) = uk37_valid.age{j}(ind);
            validData.uk37.data{i,1}(m:(m+indSum-1),1) = uk37_valid.data{j}(ind);
            validData.uk37.elev{i,1}(m:(m+indSum-1),1) = uk37_valid.elev(j);
            validData.uk37.lats{i,1}(m:(m+indSum-1),1) = uk37_valid.lats(j);
            validData.uk37.lons{i,1}(m:(m+indSum-1),1) = uk37_valid.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    allData.uk37.coreName{i,1}{k,1} = uk37_valid.coreName{j};
                    end
                else
                validData.uk37.coreName{i,1}{m:(m+indSum-1),1} = uk37_valid.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end 
    % tex86
    validData.tex86.yearRange(i,:) = binRange;
    m = 1; % initialize adjustable index
    for j = 1:length(tex86_valid.age)
        ind = tex86_valid.age{j} >= ageBin(i)+1 & tex86_valid.age{j} <= ageBin(i+1);
        indSum = nansum(ind);
        if indSum ~= 0
            % number arrays
            validData.tex86.age{i,1}(m:(m+indSum-1),1) = tex86_valid.age{j}(ind);
            validData.tex86.data{i,1}(m:(m+indSum-1),1) = tex86_valid.data{j}(ind);
            validData.tex86.elev{i,1}(m:(m+indSum-1),1) = tex86_valid.elev(j);
            validData.tex86.lats{i,1}(m:(m+indSum-1),1) = tex86_valid.lats(j);
            validData.tex86.lons{i,1}(m:(m+indSum-1),1) = tex86_valid.lons(j);
            % cell/string arrays
                if indSum > 1
                    for k = m:(m+indSum-1)
                    validData.tex86.coreName{i,1}{k,1} = tex86_valid.coreName{j};
                    end
                else
                validData.tex86.coreName{i,1}{m:(m+indSum-1),1} = tex86_valid.coreName{j};
                end
            % update m
            m = m+indSum-1;
        end
    end    
end

end
