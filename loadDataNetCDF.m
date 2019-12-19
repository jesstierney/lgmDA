function data_out = loadDataNetCDF(filename,datatype)
% function data_out = loadDataNetCDF(filename,datatype)
%
% Function to load in data from the proxy data collection. takes as input
% the filename of the netcdf that contains the data and the proxy type
% desired.
%
% Options for datatype are:
% "uk37" = UK37
% "d18o" = d18O
% "mgca" = Mg/Ca
% "tex86" = TEX86

%TESTING
filename='proxysite-2019-10-09.nc';
datatype="mgca";
%
finfo=ncinfo(filename);
%get core names
cores={finfo.Groups.Name};
%set up for uk and tex
if contains(datatype,"uk") || contains(datatype,"tex86")
    ind_1=false(1,length(cores));
else
end
%set-up for mgca and delo
if contains(datatype,"mgca") || contains(datatype,"d18o")
    ind_foram=NaN(1,length(cores));
    species_list=cell(1,length(cores));
else
end
%for mgca we also need the cleaning parameter
if contains(datatype,"mgca")
    clean_list=cell(1,length(cores));
else
end

for i=1:length(cores)
    var_list={finfo.Groups(i).Groups(2).Variables.Name};
    if contains(datatype,"uk") || contains(datatype,"tex86")
        ind_1(i)=contains(datatype,var_list);
    else
    end
    if contains(datatype,"mgca") || contains(datatype,"d18o")
        var_ind=contains(var_list,datatype);
        species_list{i}=var_list(var_ind);
        ind_foram(i)=sum(var_ind);
    else
    end
    if contains(datatype,"mgca")
        vars_now=finfo.Groups(i).Groups(2).Variables(var_ind);
        clean_now=cell(length(vars_now),1);
        for j=1:length(vars_now)
            clean_now{j}=vars_now(j).Attributes(6).Value;
        end
        clean_list{i}=clean_now;
    else
    end
end

if contains(datatype,"mgca")
    species_c=species_list(~cellfun(@isempty,species_list));
    clean_c=clean_list(~cellfun(@isempty,species_list));
elseif contains(datatype,"d18o")
    species_c=species_list(~cellfun(@isempty,species_list));
else
    species_c=cell(sum(ind_1),1);
end

%for forams, compute a logical index of where they are and record the total
if contains(datatype,"mgca") || contains(datatype,"d18o")
    ind_1=ind_foram>0;
    data_out.tot_foram=sum(ind_foram);
else
end

ages=cell(sum(ind_1),1);
proxy=cell(sum(ind_1),1);
lats=NaN(sum(ind_1),1);
lons=NaN(sum(ind_1),1);
elev=NaN(sum(ind_1),1);
core_name=cores(ind_1);
reference=cell(sum(ind_1),1);
%% get netcdf data (very slow...ncread is very slow)
tic
%read in all sites
parfor i=1:length(core_name)
    %get metadata
    refind=ismember(cores,core_name{i});
    reference{i}=finfo.Groups(refind).Attributes(7).Value;
    lats(i)=finfo.Groups(refind).Attributes(3).Value;
    lons(i)=finfo.Groups(refind).Attributes(4).Value;
    elev(i)=finfo.Groups(refind).Attributes(5).Value;
    ages{i}=ncread(filename,strcat('/',core_name{i},'/data/age_median'));
   %now get the proxy data
   if contains(datatype,"uk") || contains(datatype,"tex86")
      proxy{i}=ncread(filename,strcat('/',core_name{i},'/data/',datatype));
   elseif contains(datatype,"mgca") || contains(datatype,"d18o")
      spec_now=species_c{i};
      foram_now=cell(1,length(spec_now));
      for j=1:length(spec_now)
        foram_now{j}=ncread(filename,strcat('/',core_name{i},'/data/',spec_now{j}));
      end
      proxy{i}=foram_now;
   else
   end
end
toc
%clear finfo ans clean_now clean_list cores i ind_foram ind_1 j species_list var_ind var_list vars_now
%now put everything into a structure
data_out.ages = ages;
data_out.lats = lats;
data_out.lons = lons;
data_out.elev = elev;
data_out.reference = reference;
data_out.proxy = proxy;