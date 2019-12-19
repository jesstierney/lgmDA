%get proxies
%proxies for DA
calibfrac = .75;
%validfrac = .25;
[calibData,validData] = prepProxies (["uk","tex","mg","delo"],calibfrac,"lgm");
%[calibData,validData] = prepProxies ("uk",frac,0,"holo");

%% Initialize a cell array to hold the PSM for each site
F = cell(length(calibData.proxydata),1);
% Initialize a cell array to hold the indices for H in the PI runs
F2 = cell(length(calibData.proxydata),1);
% Initialize a cell array to for the validation points
V = cell(length(validData.proxydata),1);
% Initialize a cell array to for the validation points - PI for bias
% correction
V2 = cell(length(validData.proxydata),1);

%%
ensFile = 'holo13_50_40_cesm.ens';
biasEns = 'pi13_20_70_cesm.ens';

%ensFile = 'hol_20_75.ens';
%biasEns = 'pi_20_45.ens';

%load ensemble metadata
[ensMeta, design] = loadEnsembleMeta(ensFile);
% Load the prior ensemble which is M
M = loadEnsemble(ensFile);
%load second file as well.
M2 = loadEnsemble(biasEns);
[ensMeta2, ~] = loadEnsembleMeta(biasEns);

% get global mean SST and append to state vector
tos_meta = design.var(checkDesignVar(design, 'sst-ann'));
ocean = logical( tos_meta.meta.specs.ocean )';
ind = find(strcmp('sst-ann',ensMeta.varName));
if ind == 1
    varDex = (1:ensMeta.nEls(ind))';
else
    varDex = (sum(ensMeta.nEls(1:ind-1))+1:sum(ensMeta.nEls(1:ind)))';
end
load tarea.mat;
[nLon, nLat] = size(tarea);
tarea = reshape(tarea,nLon*nLat,1);
tarea = tarea(ocean);
globeSST = tarea'./sum(tarea) * M(varDex,:);
M = cat(1, M, globeSST );
%

% get global mean TAS and append to state vector
tas_meta = design.var(checkDesignVar(design, 'tas-ann'));
lat_cam=tas_meta.meta.lat;
lon_cam=tas_meta.meta.lon;
ind = find(strcmp('tas-ann',ensMeta.varName));
if ind == 1
    varDex = (1:ensMeta.nEls(ind))';
else
    varDex = (sum(ensMeta.nEls(1:ind-1))+1:sum(ensMeta.nEls(1:ind)))';
end
gw = repmat(cos(lat_cam * pi / 180),1,144)';
gw = reshape(gw,length(lon_cam)*length(lat_cam),1);
globeTAS = gw' ./ sum(gw) * M(varDex,:);
M = cat(1, M, globeTAS );
%
monthName = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
% set bayes files
filesdelo=["poolann_params.mat";"poolsea_params.mat";"hiersea_params_cut.mat"];
filesmg = ["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params_cut.mat"];
fileuk = 'bayes_posterior_clip.mat';
%% Initialize the PSM for each site.
tic
parfor s = 1:length(calibData.proxydata)
    if strcmp(calibData.proxytype{s},'uk')
        F{s} = ukPSM(calibData.proxylats(s),calibData.proxylons(s),'bayesFile',fileuk);
        F{s}.getStateIndices( ensMeta, "sst-mon", monthName );
        [Xt,~,~] = modernTSd18O(calibData.proxylats(s),calibData.proxylons(s));
        Xt = Xt(:);
        F2{s}.H = getClosestLatLonIndex([calibData.proxylats(s) calibData.proxylons(s)],ensMeta2,"sst-mon", "time", monthName);
    elseif strcmp(calibData.proxytype{s},'tex')
        F{s} = texPSM(calibData.proxylats(s),calibData.proxylons(s));
        F{s}.getStateIndices( ensMeta, "sst-ann");
        [Xt,~,~] = modernTSd18O(calibData.proxylats(s),calibData.proxylons(s));
        Xt = mean(Xt);
        F2{s}.H = getClosestLatLonIndex([calibData.proxylats(s) calibData.proxylons(s)],ensMeta2,"sst-ann");
    elseif strcmp(calibData.proxytype{s},'mg')
        F{s} = mgcaPSM(calibData.proxylats(s),calibData.proxylons(s),calibData.ph(s),calibData.omega(s),calibData.clean(s),calibData.species(s),'Bayes',filesmg);
        F{s}.getStateIndices( ensMeta, "sst-mon","sss-mon", monthName);
        [tempt,temps,~] = modernTSd18O(calibData.proxylats(s),calibData.proxylons(s));
        Xt = [tempt'; temps'];
        F2{s}.H = getClosestLatLonIndex([calibData.proxylats(s) calibData.proxylons(s)],ensMeta2,["sst-mon";"sss-mon"], "time", monthName);
    elseif strcmp(calibData.proxytype{s},'delo')
        F{s} = deloPSM(calibData.proxylats(s),calibData.proxylons(s),calibData.species(s),'Bayes',filesdelo);
        F{s}.getStateIndices( ensMeta, "sst-mon","d18osw-ann", monthName);
        [tempt,~,tempd] = modernTSd18O(calibData.proxylats(s),calibData.proxylons(s));
        Xt = [tempt'; tempd'];
        F2{s}.H = getClosestLatLonIndex([calibData.proxylats(s) calibData.proxylons(s)],ensMeta2,"sst-mon", "time", monthName);
        F2{s}.H(13) = getClosestLatLonIndex([calibData.proxylats(s) calibData.proxylons(s)], ensMeta2, "d18osw-ann");
    else
        error('Proxy type not recognized')
    end
    %get the constant correction
    %addConstant = getMeanAdjustment( Xt,  M2(F2{s}.H,:));
    % Activate the bias corrector
    %F{s}.useMeanCorrector( addConstant );
end
toc

%% set up Posterior Ye for validation sites
tic
parfor s = 1:length(validData.proxydata)
    if strcmp(validData.proxytype{s},'uk')
        V{s} = ukPSM(validData.proxylats(s),validData.proxylons(s));
        V{s}.getStateIndices( ensMeta, "sst-mon", monthName );
        [Xt,~,~] = modernTSd18O(validData.proxylats(s),validData.proxylons(s));
        Xt = Xt(:);
        V2{s}.H = getClosestLatLonIndex([validData.proxylats(s) validData.proxylons(s)],ensMeta2,"sst-mon", "time", monthName);
    elseif strcmp(validData.proxytype{s},'tex')
        V{s} = texPSM(validData.proxylats(s),validData.proxylons(s));
        V{s}.getStateIndices( ensMeta, "sst-ann");
        [Xt,~,~] = modernTSd18O(validData.proxylats(s),validData.proxylons(s));
        Xt = mean(Xt);
        V2{s}.H = getClosestLatLonIndex([validData.proxylats(s) validData.proxylons(s)],ensMeta2,"sst-ann");
    elseif strcmp(validData.proxytype{s},'mg')
        V{s} = mgcaPSM(validData.proxylats(s),validData.proxylons(s),validData.ph(s),validData.omega(s),validData.clean(s),validData.species(s));
        V{s}.getStateIndices( ensMeta, "sst-mon","sss-mon", monthName);
        [tempt,temps,~] = modernTSd18O(validData.proxylats(s),validData.proxylons(s));
        Xt = [tempt'; temps'];
        V2{s}.H = getClosestLatLonIndex([validData.proxylats(s) validData.proxylons(s)],ensMeta2,["sst-mon";"sss-mon"], "time", monthName);
    elseif strcmp(validData.proxytype{s},'delo')
        V{s} = deloPSM(validData.proxylats(s),validData.proxylons(s),validData.species(s));
        V{s}.getStateIndices( ensMeta, "sst-mon","d18osw-ann", monthName);
        [tempt,~,tempd] = modernTSd18O(validData.proxylats(s),validData.proxylons(s));
        Xt = [tempt'; tempd'];
        V2{s}.H = getClosestLatLonIndex([validData.proxylats(s) validData.proxylons(s)],ensMeta2,"sst-mon", "time", monthName);
        V2{s}.H(13) = getClosestLatLonIndex([validData.proxylats(s) validData.proxylons(s)], ensMeta2, "d18osw-ann");
    else
        error('Proxy type not recognized')
    end
    %get the constant correction
    %addConstant = getMeanAdjustment( Xt,  M2(V2{s}.H,:));
    % Activate the bias corrector
    %V{s}.useMeanCorrector( addConstant );
end
toc
%% run the assimilation!
% set D the proxy data
D = calibData.proxydata;
% set R the proxy variance. If you set this to NaN then it comes from the forward model.
R = NaN(length(calibData.proxydata),1);
% other option is to specify R specifically:
factorn = 100;
R(calibData.proxytype == "uk") = (0.05^2)/factorn; %mean posterior spread for UK model
R(calibData.proxytype == "tex") = (0.05^2)/factorn; %typical posterior spread for TEX86 model
R(calibData.proxytype == "mg") = (0.13^2)/factorn; %mean RMSE from BayMAG in log units
R(calibData.proxytype == "delo") = (0.5^2)/factorn;  %mean posterior spread for d18O model
%localization
[w, yloc] = covLocalization([calibData.proxylats calibData.proxylons], ensMeta, 12000 );
%need to add on rows to w
w(end+1,:)=1;
w(end+1,:)=1;
%%
tic
[Amean, Avar, Ye, Rpost, ~, calib] = dash( M, D, R, F,'localize',{w, yloc});
toc
%% regrid
gridSize = tos_meta.meta.specs.gridSize;
sst_ann = regridTripolar( Amean, 'sst-ann', ensMeta, design, ocean, gridSize);
sst_ann_var = regridTripolar( Avar, 'sst-ann', ensMeta, design, ocean, gridSize);
prior_sst_ann = regridTripolar(M, 'sst-ann', ensMeta, design, ocean, gridSize);
tas_ann = regridAnalysis( Amean, 'tas-ann', ensMeta, design);
pr_ann = regridAnalysis( Amean, 'pr-ann', ensMeta, design);
prior_tas_ann = regridAnalysis(M, 'tas-ann', ensMeta, design);

pop_locs=reshape(tos_meta.meta.tri,size(sst_ann,1),size(sst_ann,2),2);
lat_pop = pop_locs(:,:,1);
lon_pop = pop_locs(:,:,2);

GSST = Amean(size(Amean,1)-1,:);
GAT = Amean(size(Amean,1),:)-273.15;
varDex = (1:ensMeta.nEls(1))';
globeSSTpost = tarea'./sum(tarea) * Amean(varDex);
ind = find(strcmp('tas-ann',ensMeta.varName));
varDex = (sum(ensMeta.nEls(1:ind-1))+1:sum(ensMeta.nEls(1:ind)))';
globeTASpost = gw' ./ sum(gw) * Amean(varDex) - 273.15;
%% validation
tot = Amean + Avar;
Ye_valid = NaN(length(validData.proxydata),size(M,2));
Ye_prior = NaN(length(validData.proxydata),size(M,2));
R_valid = NaN(length(validData.proxydata),1);
for s=1:length(validData.proxydata)
    % Get the model values being passed
    Mpsm = tot(V{s}.H,:);
    Mp = M(V{s}.H,:);
    %Mpsm = Amean(V{s}.H) + V{s}.bias.addConstant;
    [Ye_valid(s,:), R_valid(s)] = V{s}.runForwardModel( Mpsm, NaN, 1 );
    [Ye_prior(s,:), ~] = V{s}.runForwardModel( Mp, NaN, 1 );
end

%normalize data by proxy type
proxy_norm=NaN(length(validData.proxydata),1);
Ye_norm=NaN(length(validData.proxydata),size(M,2));
Ye_prior_norm=NaN(length(validData.proxydata),size(M,2));
uk_ind=contains(validData.proxytype,"uk");
mg_ind=contains(validData.proxytype,"mg");
delo_ind=contains(validData.proxytype,"delo");
tex_ind=contains(validData.proxytype,"tex");
%
proxy_norm(uk_ind)=normalize(validData.proxydata(uk_ind));
proxy_norm(mg_ind)=normalize(validData.proxydata(mg_ind));
proxy_norm(delo_ind)=-normalize(validData.proxydata(delo_ind));
proxy_norm(tex_ind)=normalize(validData.proxydata(tex_ind));
%
Ye_norm(uk_ind,:)=(Ye_valid(uk_ind,:) - mean(validData.proxydata(uk_ind)))./std(validData.proxydata(uk_ind));
Ye_norm(mg_ind,:)=(Ye_valid(mg_ind,:) - mean(validData.proxydata(mg_ind)))./std(validData.proxydata(mg_ind));
Ye_norm(delo_ind,:)=-(Ye_valid(delo_ind,:) - mean(validData.proxydata(delo_ind)))./std(validData.proxydata(delo_ind));
Ye_norm(tex_ind,:)=(Ye_valid(tex_ind,:) - mean(validData.proxydata(tex_ind)))./std(validData.proxydata(tex_ind));
%
Ye_prior_norm(uk_ind,:)=(Ye_prior(uk_ind,:) - mean(validData.proxydata(uk_ind)))./std(validData.proxydata(uk_ind));
Ye_prior_norm(mg_ind,:)=(Ye_prior(mg_ind,:) - mean(validData.proxydata(mg_ind)))./std(validData.proxydata(mg_ind));
Ye_prior_norm(delo_ind,:)=-(Ye_prior(delo_ind,:) - mean(validData.proxydata(delo_ind)))./std(validData.proxydata(delo_ind));
Ye_prior_norm(tex_ind,:)=(Ye_prior(tex_ind,:) - mean(validData.proxydata(tex_ind)))./std(validData.proxydata(tex_ind));
%stats
validR2 = rsquare(proxy_norm,mean(Ye_norm,2));
validRMSE = rmse(proxy_norm,mean(Ye_norm,2));
validCE = mean(1-sum((proxy_norm - Ye_norm).^2)./sum((proxy_norm - mean(proxy_norm)).^2));
calibCR = (mean(Ye,2) - calibData.proxydata).^2 ./ (var(Ye,[],2) + Rpost);
validCR = (mean(Ye_valid,2) - validData.proxydata).^2 ./ (var(Ye_valid,[],2) + R_valid);
varRatio = Rpost./var(Ye,0,2);
validR = corr(proxy_norm,mean(Ye_norm,2));
%quick plot
f2=figure(2); clf;
set(f2,'pos',[600 500 400 400]);
p1=scatter(proxy_norm(uk_ind),mean(Ye_norm(uk_ind,:),2)); hold on;
p2=scatter(proxy_norm(mg_ind),mean(Ye_norm(mg_ind,:),2));
p3=scatter(proxy_norm(delo_ind),mean(Ye_norm(delo_ind,:),2));
p4=scatter(proxy_norm(tex_ind),mean(Ye_norm(tex_ind,:),2));
line([-2.5 2.5],[-2.5 2.5],'color','k');
set(gca,'xlim',[-2.5 2.5],'ylim',[-2.5 2.5]);
legend([p1 p2 p3 p4],'uk','mg','delo','tex','location','southeast');
xlabel('Proxy value');
ylabel('Validation Ye');

%%


%make Avar an average variance
AvarAvg = mean(Avar.^2,2);
runnotes = {'40 50-yr ensemble members, no inflation, no bias correction, localization 12000, R /100'};
%% save stuff
save('holo13_da_results_50_1_R100.mat','Amean','AvarAvg','calibCR','validCR','Ye','Rpost','varRatio','calibData','validR2','validCE','validData','Ye_valid','Ye_norm','proxy_norm','runnotes');