function [calibData, validData] = prepProxies (types,fCalib,timeInt)
%types a string array, e.g. types = ["uk","mg","delo"];
%fCalib and fValid are fractions of total data to use. 
%timeInt = "holo" or "lgm"
listTypes=["uk","mg","delo","tex"];
listTimes=["holo","lgm"];
%nProxCalib = round(fCalib/length(types));
%nProxValid = round(fValid/length(types));

if any (~ismember(types,listTypes))
    error('Proxy type not recognized');
else
end

if ~ismember(timeInt,listTimes)
    error('Time interval not recognized');
else
end

ct=0;

for i=1:length(types)
    if strcmp("uk",types(i))
        filename="uks_forda.mat";
        if strcmp(timeInt,"holo")
            load(filename,'lat_h','lon_h','uk_h');
            proxylatsn = lat_h;
            proxylonsn = lon_h;
            proxydatan = uk_h;
        else
            load(filename,'lats_lgm','lons_lgm','uk_lgm');
            proxylatsn = lats_lgm;
            proxylonsn = lons_lgm;
            proxydatan = uk_lgm;
        end
    elseif strcmp("tex",types(i))
        filename="tex_forda.mat";
        if strcmp(timeInt,"holo")
            load(filename,'lat_h','lon_h','tex_holo');
            proxylatsn = lat_h;
            proxylonsn = lon_h;
            proxydatan = tex_holo;
        else
            load(filename,'lats_lgm','lons_lgm','tex_lgm');
            proxylatsn = lats_lgm;
            proxylonsn = lons_lgm;
            proxydatan = tex_lgm;
        end
    elseif strcmp("mg",types(i))
        filename="mg_forda.mat";
        if strcmp(timeInt,"holo")
            load(filename,'lat_h','lon_h','mg_h','clean_h','omega_h','ph_h','species_h');
            proxylatsn = lat_h;
            proxylonsn = lon_h;
            proxydatan = log(mg_h);
            cleann = clean_h;
            omegan = omega_h;
            phn = ph_h;
            speciesn = species_h;
        else
            load(filename,'lats_lgm','lons_lgm','mg_lgm','clean_lgm','omega_lgm','ph_lgm_adj','species_lgm');
            proxylatsn = lats_lgm;
            proxylonsn = lons_lgm;
            proxydatan = log(mg_lgm);
            cleann = clean_lgm;
            omegan = omega_lgm;
            phn = ph_lgm_adj;
            speciesn = species_lgm;    
        end
    else %type is delo
        filename="delo_forda.mat";
        if strcmp(timeInt,"holo")
            load(filename,'lat_h','lon_h','delo_h','species_h','dsw_h');
            proxylatsn = lat_h;
            proxylonsn = lon_h;
            proxydatan = delo_h;
            speciesn = species_h;
            dswn = dsw_h;
        else
            load(filename,'lats_lgm','lons_lgm','delo_lgm','species_lgm','dsw_lgm');
            proxylatsn = lats_lgm;
            proxylonsn = lons_lgm;
            proxydatan = delo_lgm;
            speciesn = species_lgm;
            dswn = dsw_lgm;
        end
    end
    %downsample
    %common to all
    Np = length(proxylatsn);
    nProxCalib = round(Np * fCalib);
    nProxValid = Np - nProxCalib;
    inds = randsample(Np, nProxCalib + nProxValid);
    proxylatsc = proxylatsn(inds(1:nProxCalib));
    proxylonsc = proxylonsn(inds(1:nProxCalib));
    proxydatac = proxydatan(inds(1:nProxCalib));
    proxylatsv = proxylatsn(inds(nProxCalib+1:end));
    proxylonsv = proxylonsn(inds(nProxCalib+1:end));
    proxydatav = proxydatan(inds(nProxCalib+1:end));
    proxytypec=repmat(types(i),length(proxydatac),1);
    proxytypev=repmat(types(i),length(proxydatav),1);
    %append
    if i==1
        calibData.proxylats = proxylatsc;
        calibData.proxylons = proxylonsc;
        calibData.proxydata = proxydatac;
        calibData.proxytype = proxytypec;
        validData.proxylats = proxylatsv;
        validData.proxylons = proxylonsv;
        validData.proxydata = proxydatav;
        validData.proxytype = proxytypev;
    else
        calibData.proxylats = [calibData.proxylats; proxylatsc];
        calibData.proxylons = [calibData.proxylons; proxylonsc];
        calibData.proxydata = [calibData.proxydata; proxydatac];
        calibData.proxytype = [calibData.proxytype; proxytypec];
        validData.proxylats = [validData.proxylats; proxylatsv];
        validData.proxylons = [validData.proxylons; proxylonsv];
        validData.proxydata = [validData.proxydata; proxydatav];
        validData.proxytype = [validData.proxytype; proxytypev];
    end
    forams = strcmp(types(i),"mg") || strcmp(types(i),"delo");
    if forams
        ct = ct+1;
        speciesc = speciesn(inds(1:nProxCalib));
        speciesv = speciesn(inds(nProxCalib+1:end));
        if ct <= 1
        speciestc = speciesc;
        speciestv = speciesv;
        else
        speciestc = [speciestc; speciesc];
        speciestv = [speciestv; speciesv];
        end
    else
    end
    if strcmp("mg",types(i))
        cleanc = cleann(inds(1:nProxCalib));
        omegac = omegan(inds(1:nProxCalib));
        phc = phn(inds(1:nProxCalib));
        cleanv = cleann(inds(nProxCalib+1:end));
        omegav = omegan(inds(nProxCalib+1:end));
        phv = phn(inds(nProxCalib+1:end));
    else
    end
    if strcmp("delo",types(i))
        dswc = dswn(inds(1:nProxCalib));
        dswv = dswn(inds(nProxCalib+1:end));
    else
    end
end
if forams
    calibData.species = strings(length(calibData.proxydata),1);
    calibData.species(strcmp("mg",calibData.proxytype) | strcmp("delo",calibData.proxytype)) = speciestc;
    validData.species = strings(length(validData.proxydata),1);
    validData.species(strcmp("mg",validData.proxytype) | strcmp("delo",validData.proxytype)) = speciestv;
else
end
%mg is the only proxy with clean, omega, ph
if ismember("mg",types)
    calibData.omega = NaN(length(calibData.proxydata),1);
    calibData.ph = NaN(length(calibData.proxydata),1);
    calibData.clean = NaN(length(calibData.proxydata),1);
    validData.omega = NaN(length(validData.proxydata),1);
    validData.ph = NaN(length(validData.proxydata),1);
    validData.clean = NaN(length(validData.proxydata),1);

    calibData.omega(strcmp("mg",calibData.proxytype))=omegac;
    calibData.ph(strcmp("mg",calibData.proxytype))=phc;
    calibData.clean(strcmp("mg",calibData.proxytype))=cleanc;
    validData.omega(strcmp("mg",validData.proxytype))=omegav;
    validData.ph(strcmp("mg",validData.proxytype))=phv;
    validData.clean(strcmp("mg",validData.proxytype))=cleanv;
else
end

%delo is the only proxy with delosw, needed for models w/o this variable.
if ismember("delo",types)
    calibData.dsw = NaN(length(calibData.proxydata),1);
    validData.dsw = NaN(length(validData.proxydata),1);

    calibData.dsw(strcmp("delo",calibData.proxytype))=dswc;
    validData.dsw(strcmp("delo",validData.proxytype))=dswv;
else
end

%sanity check plot
% figure(2); clf;
% m_proj('mollweide','lon',[-180 180]); hold on;
% m_coast('patch',[.6 .6 .6]);
% s1=m_scatter(calibData.proxylons,calibData.proxylats,50,'filled');
% s2=m_scatter(validData.proxylons,validData.proxylats,50,'filled');
% m_grid('xtick',[],'ytick',[]);
% legend([s1 s2],'Calibration','Validation');