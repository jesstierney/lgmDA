clear
home_fold = cd;
cd ../; addpath(genpath('DASH/')); cd DASH/

processNetCDF = false; % warning! Should only need to run this once; If this is set to true, be prepared to wait forevaaaa

% Written by M. Osman, Jan2020 (mattosman@arizona.edu) in Matlab2018b via DASH v3.2.x
% This script maninpulates five iCESM *.nc fields (see below; must state "processNetCDF = true, above)
% into a format that allows DASH to create .grid mean (decadal-)monthly and mean (decadal-) annual instructions. 
% Accomponying this 
%   To run, requires the (iCESM) files: 
%     1. d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc *
%     2. sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc *
%     3. tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc *
%     4. toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc *
%     5. tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc *
%       NOTE* If using on a different desktop, you will need to route the 
%       script to these files by changing the cd() instructions in lines 89 + 131 !
%     6. A subfolder in your Current Folder called ModelData/for outputting .mat files 
%   
%% Several existing data-wrangling options within the Tierney Lab, identified till date (Jan 2020)
%     % iCESM decadal-average (for each month) time slice simulations compiled by Brewster Malevich - these configs are useful! 
%     cd('/Volumes/Osman_climate_data_repository/LMR/data/model/icesm_lgm2present');
%     cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present'); % or, (THIS IS WHAT THIS CODE BASE USES!) 
% 
%     % can also get the original iCESM data from here: note, 12ky is missing(?) 
%     cd('/Volumes/Lab/Projects/DTDA/icesm/'); 
% 
%     % TRACE21k simulation compiled by Brewster Malevich - these are decadal averages only!
%     cd('/Volumes/Osman_climate_data_repository/LMR/data/model/ccsm3_trace21ka');
% 
%     % TRACE21k pop (ocean model) monthly; connected to Lab server (only sst and sss)
%     % the file of interest is the big one (trace_pop_temp_400cm.nc) - use
%     %   formatting of read_trace.m to break into small pieces, if necessary
%     cd('/Volumes/Lab/Projects/DTDA/trace_monthly/');
% 
%     % TRACE21k cam (atm model) monthly; also,connected to Lab server; (surf air temp)
%     % (decadal pop averages of SST and SSS also available in this folder, oddly) 
%     cd('/Volumes/Lab/Projects/DTDA/trace/');

%% Check if .grid files already exist!  If so, throw up an option to overwrite and continue

if exist('iCESM_MonAvg_Ocean.grid') > 0 
    answer = questdlg('The file ''iCESM_MonAvg_Ocean.grid'' exist in your currently folder! Overwrite it?','Overwrite file?','Yes','No','No');
    switch answer
        case 'Yes'
            deleteOpt = true;
        case 'No'
            deleteOpt = false; end
    if deleteOpt
    recycle('on'); delete('iCESM_MonAvg_Ocean.grid'); recycle('off'); else; return; end
elseif exist('iCESM_AnnAvg_Ocean.grid') > 0 
    answer = questdlg('The file ''iCESM_AnnAvg_Ocean.grid'' exist in your currently folder! Overwrite it?','Overwrite file?','Yes','No','No');
    switch answer
        case 'Yes'
            deleteOpt = true;
        case 'No'
            deleteOpt = false; end
    if deleteOpt
    recycle('on'); delete('iCESM_AnnAvg_Ocean.grid'); recycle('off'); else; return; end
elseif exist('iCESM_MonAvg_Atm.grid') > 0 
    answer = questdlg('The file ''iCESM_MonAvg_Atm.grid'' exist in your currently folder! Overwrite it?','Overwrite file?','Yes','No','No');
    switch answer
        case 'Yes'
            deleteOpt = true;
        case 'No'
            deleteOpt = false; end
    if deleteOpt
    recycle('on'); delete('iCESM_MonAvg_Atm.grid'); recycle('off'); else; return; end
elseif exist('iCESM_AnnAvg_Atm.grid') > 0 
    answer = questdlg('The file ''iCESM_AnnAvg_Atm.grid'' exist in your currently folder! Overwrite it?','Overwrite file?','Yes','No','No');
    switch answer
        case 'Yes'
            deleteOpt = true;
        case 'No'
            deleteOpt = false; end
    if deleteOpt
    recycle('on'); delete('iCESM_AnnAvg_Atm.grid'); recycle('off'); else; return; end
end

%% Need to first convert the iCESM brewster .nc files to .mat files
% % why? because Jon's DASH programming doesn't support a month dimension.
% % Note I could alternatively reshape the matrix to be 1 less dimension by
% % putting the months - years as 1 single dimension.  

if processNetCDF
    
    % ANNUAL 
    cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present');

            slice = ["iPI.01";"i03ka.01";"i06ka.03";"i09ka.01";"i12ka.01";"i18ka.01";"i21ka.03"];    

            lat = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            d18Osw = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','d18Osw');
                d18Osw = squeeze(nanmean(d18Osw,3));
            save d18Osw.mat d18Osw lat lon time slice -v7.3

            lat = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            sos = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','sos');
                sos = squeeze(nanmean(sos,3));
            save sos.mat sos lat lon time slice -v7.3

            lat = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            tas = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','tas');
                tas = squeeze(nanmean(tas,3));
            save tas.mat tas lat lon time slice -v7.3   

            lat = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            toGA = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','toGA');
                toGA = squeeze(nanmean(toGA,3));
            save toGA.mat toGA lat lon time slice -v7.3   

            lat = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            tos = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','tos');
                tos = squeeze(nanmean(tos,3));
            save tos.mat tos lat lon time slice -v7.3   

    cd(home_fold);

    % MONTHLY 
    cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present');

            slice = ["iPI.01";"i03ka.01";"i06ka.03";"i09ka.01";"i12ka.01";"i18ka.01";"i21ka.03"];    

            lat = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            d18Osw = ncread('d18Osw_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','d18Osw');
            d18Osw = reshape(d18Osw,size(d18Osw,1),size(d18Osw,2),size(d18Osw,3)*size(d18Osw,4),size(d18Osw,5)); 
            time_month = ((0.5:11.5)./12 + repmat(time-1,1,12))'; time = time_month(:);
            save d18Osw_monthly.mat d18Osw lat lon time slice -v7.3

            lat = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            sos = ncread('sos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','sos');
            sos = reshape(sos,size(sos,1),size(sos,2),size(sos,3)*size(sos,4),size(sos,5)); 
            time_month = ((0.5:11.5)./12 + repmat(time-1,1,12))'; time = time_month(:);
            save sos_monthly.mat sos lat lon time slice -v7.3

            lat = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            tas = ncread('tas_sfc_Adecmon_iCESM_LGMtoPresent_000101-090012.nc','tas');
            tas = reshape(tas,size(tas,1),size(tas,2),size(tas,3)*size(tas,4),size(tas,5)); 
            time_month = ((0.5:11.5)./12 + repmat(time-1,1,12))'; time = time_month(:);
            save tas_monthly.mat tas lat lon time slice -v7.3

            lat = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            toGA = ncread('toGA_0-200m_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','toGA');
            toGA = reshape(toGA,size(toGA,1),size(toGA,2),size(toGA,3)*size(toGA,4),size(toGA,5)); 
            time_month = ((0.5:11.5)./12 + repmat(time-1,1,12))'; time = time_month(:);
            save toGA_monthly.mat toGA lat lon time slice -v7.3

            lat = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lat'); 
            lon = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','lon');
            time = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','time')./365; 
            tos = ncread('tos_sfc_Odecmon_iCESM_LGMtoPresent_000101-090012.nc','tos');
            tos = reshape(tos,size(tos,1),size(tos,2),size(tos,3)*size(tos,4),size(tos,5)); 
            time_month = ((0.5:11.5)./12 + repmat(time-1,1,12))'; time = time_month(:);
            save tos_monthly.mat tos lat lon time slice -v7.3

    cd(home_fold);

end

%% Step1: create .grid file 
% must create two separate .grid files, to account for different gridding

if createGridFile
% REMINDER: possible .grid dimensions:
%   lon: Longitude / x-axis
%   lat: Latitude / y-axis
%   lev: Height / z-axis
%   tri: A tripolar dimension (More on this later in the tutorial.)
%   time: A time dimension
%   run: A model-ensemble dimension
%   var: A dimension for different variables.

% Annual 

    % ocean (need to split up from atmosphere because tripolar gridding
    cd ModelData/
        load d18Osw.mat % contains d18Osw - lat - lon - slice - time (decadal monthly averages)
        load sos.mat    % contains sos - lat - lon - slice - time (decadal monthly averages)
        load tos.mat    % contains tos - lat - lon - slice - time (decadal monthly averages)
    cd ../
        tri = [lat(:),lon(:)]; % create tripolar gridding
        % now create grid file
        variables = ["d18Osw";"sos";"tos"];
        dimOrder = ["tri","tri","time","run"];
    % create new .grid file
        meta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', variables );
        grid = gridFile.new('iCESM_AnnAvg_Ocean.grid', meta );
    % now add .mat instructions
    cd ModelData/
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'd18Osw' ); % redefine meta before adding d18Osw
        grid.addData( 'mat', 'd18Osw.mat', 'd18Osw', dimOrder, sourceMeta );
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'sos' ); % redefine meta before adding sos
        grid.addData( 'mat', 'sos.mat', 'sos', dimOrder, sourceMeta );
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'tos' ); % redefine meta before adding tos
        grid.addData( 'mat', 'tos.mat', 'tos', dimOrder, sourceMeta );
    cd ../

    % atmosphere 
    cd ModelData/
        load tas.mat % contains tas - lat - lon - slice - time (decadal monthly averages); 
    cd ../
        variables = ["tas"];
        dimOrder = ["lon","lat","time","run"];
    % create new grid file
        meta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'time', time, 'run', slice,'var', variables );
        grid = gridFile.new('iCESM_AnnAvg_Atm.grid', meta );
    % now add .mat instructions
    cd ModelData/
        sourceMeta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'time', time, 'run', slice,'var', 'tas' ); % redefine meta before adding d18Osw
        grid.addData( 'mat', 'tas.mat', 'tas', dimOrder, sourceMeta );
    cd ../

% Monthly

    % ocean (need to split up from atmosphere because tripolar gridding
    cd ModelData/
        load d18Osw_monthly.mat % contains d18Osw - lat - lon - slice - time (decadal monthly averages)
        load sos_monthly.mat    % contains sos - lat - lon - slice - time (decadal monthly averages)
        load tos_monthly.mat    % contains tos - lat - lon - slice - time (decadal monthly averages)
    cd ../
        tri = [lat(:),lon(:)]; % create tripolar gridding
        % now create grid file
        variables = ["d18Osw";"sos";"tos"];
        dimOrder = ["tri","tri","time","run"];
    % create new .grid file
        meta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', variables );
        grid = gridFile.new('iCESM_MonAvg_Ocean.grid', meta );
    % now add .mat instructions
    cd ModelData/
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'd18Osw' ); % redefine meta before adding d18Osw
        grid.addData( 'mat', 'd18Osw_monthly.mat', 'd18Osw', dimOrder, sourceMeta );
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'sos' ); % redefine meta before adding sos
        grid.addData( 'mat', 'sos_monthly.mat', 'sos', dimOrder, sourceMeta );
        sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice,'var', 'tos' ); % redefine meta before adding tos
        grid.addData( 'mat', 'tos_monthly.mat', 'tos', dimOrder, sourceMeta );
    cd ../

    % atmosphere 
    cd ModelData/
        load tas_monthly.mat % contains tas - lat - lon - slice - time (decadal monthly averages); 
    cd ../
        variables = ["tas"];
        dimOrder = ["lon","lat","time","run"];
    % create new grid file
        meta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'time', time, 'run', slice,'var', variables );
        grid = gridFile.new('iCESM_MonAvg_Atm.grid', meta );
    % now add .mat instructions
    cd ModelData/
        sourceMeta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'time', time, 'run', slice,'var', 'tas' ); % redefine meta before adding d18Osw
        grid.addData( 'mat', 'tas_monthly.mat', 'tas', dimOrder, sourceMeta );
    cd ../
    
end

%% Step1: create .grid file
% % this is option 2, combinging atmosphere and ocean into one grid file
% % NOTE - THIS DOESN'T APPEAR TO WORK!  Can't combine different file sizes into a single .grid file!
%   
% % define all variables fist
% variables = ["d18Osw";"sos";"tos";"tas"];
% 
% % load ocean var's
% cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present');
%     load d18Osw % contains sos lat lon slice time (decadal averages)
%     load sos % contains sos lat lon slice time (decadal averages)
%     load tos % contains sos lat lon slice time (decadal averages)
%     tri = [lat(:),lon(:)]; % create tripolar gridding (need to do this here because following line loads over w/ non-tri lat and lon data) 
%     load tas % contains tas lat lon slice time (decadal averages)
% cd(home_fold);    
%     
% % add ocean directions to .grid
%     dimOrder = ["tri","tri","time","run"];
%     % now create grid file
%     meta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'tri', tri, 'time', time, 'run', slice,'var', variables );
%     grid = gridFile.new('iCESM_DecAvg.grid', meta );
%     % append ocean .mat files
%     cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present');
%         sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice, 'var', 'd18Osw' ); % redefine meta before adding d18Osw
%         grid.addData( 'mat', 'd18Osw.mat', 'd18Osw', dimOrder, sourceMeta );
%         sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice, 'var', 'sos' ); % redefine meta before adding sos
%         grid.addData( 'mat', 'sos.mat', 'sos', dimOrder, sourceMeta );
%         sourceMeta = gridFile.defineMetadata('tri', tri, 'time', time, 'run', slice, 'var', 'tos' ); % redefine meta before adding tos
%         grid.addData( 'mat', 'tos.mat', 'tos', dimOrder, sourceMeta );
%     cd(home_fold);
% 
% % add atm. directions to .grid
%     dimOrder = ["lon","lat","time","run"];
%     cd('/Users/matthewosman/LMR/data/model/icesm_lgm2present');
%         sourceMeta = gridFile.defineMetadata('lon', lon, 'lat', lat, 'time', time, 'run', slice, 'var', 'tas' ); % redefine meta before adding d18Osw
%         grid.addData( 'mat', 'tas.mat', 'tas', dimOrder, sourceMeta );
%     cd(home_fold);
