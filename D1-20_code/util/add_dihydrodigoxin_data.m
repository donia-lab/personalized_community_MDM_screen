function manifest = add_dihydrodigoxin_data(manifest,dihydrodigoxin_loc,DMSO)

%This function takes an already populated manifest and the data from manual
%dihydrodigoxin data and adds the manual data to the manifest

data = readtable(dihydrodigoxin_loc,'ReadRowNames',false,...
    'ReadVariableNames',true);
data.name = matlab.lang.makeValidName(data.name);

[~,idata,imanifest] = intersect(data.name,manifest.name);

for i = 1:length(idata)
    met_level = data.digoxin_met_1(idata(i));
    
    if DMSO
        manifest.digoxin_met_1(imanifest(i)) = met_level;
        manifest.norm_digoxin_met_1(imanifest(i)) = ...
            met_level/manifest.istd(imanifest(i));
    else
        manifest.met_1(imanifest(i)) = met_level;
        manifest.norm_met_1(imanifest(i)) = ...
            met_level/manifest.istd(imanifest(i));
    end
    
end

