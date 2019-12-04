function mod_manifest = import_16S_results(asv_file,species_file,...
    genus_file,family_file,order_file,manifest_file,primary_biomass_file,secondary_biomass_file)
%This function takes a manifest and 16S result file locations and returns a
%manifest with the community compositions as table entries

%Turn off warning for valid name modification, the script deals with this
%internally by using matlab.lang.makeValidName
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%Import tables and manifests
asv_table = readtable(asv_file,'ReadRowNames',true);
order_table = readtable(order_file,'ReadRowNames',true);
family_table = readtable(family_file,'ReadRowNames',true);
genus_table = readtable(genus_file,'ReadRowNames',true);
species_table = readtable(species_file,'ReadRowNames',true);
manifest = readtable(manifest_file,'ReadRowNames',false);

mod_manifest = manifest; 

sample_names = manifest.sample;

num_sample = length(sample_names); 

asv_cell = cell(num_sample,1);
rel_asv_cell = asv_cell;
order_cell = asv_cell;
family_cell = asv_cell;
genus_cell = asv_cell;
species_cell = asv_cell;

for i = 1:num_sample
    
    %Get matlab valid name to control for matlab's import behavior
    orig_name = sample_names{i};
    valid_name = matlab.lang.makeValidName(orig_name);
    
    %Add entries to cells from the raw tables
    asv_cell{i} = asv_table(:,valid_name); 
    rel_asv_table = asv_cell{i};
    rel_asv_table.Variables = rel_asv_table.Variables/sum(rel_asv_table.Variables); 
    rel_asv_cell{i} = rel_asv_table;
    
    order_cell{i} = order_table(:,valid_name);
    family_cell{i} = family_table(:,valid_name); 
    genus_cell{i} = genus_table(:,valid_name);
    species_cell{i} = species_table(:,valid_name);
end

%Add taxonomical data to the manifest
mod_manifest.asv = asv_cell; 
mod_manifest.rel_asv = rel_asv_cell;
mod_manifest.order = order_cell;
mod_manifest.family = family_cell;
mod_manifest.genus = genus_cell;
mod_manifest.species = species_cell;

%Import the biomass into the manifest
media_list = unique(mod_manifest.media);
media_list = media_list(~contains(media_list,{'feces','none'}));

biomass_data = readtable(primary_biomass_file);
secondary_biomass_data = readtable(secondary_biomass_file);

for i = 1:size(mod_manifest,1)
   sample_media = mod_manifest.media{i};
   sample_donor = mod_manifest.donor(i);
   if contains(sample_media,media_list)

      matching_mass = (contains(biomass_data.media,sample_media)) & ...
          (biomass_data.donor == sample_donor);
            
      mod_manifest.biomass(i) =mean(biomass_data.wet_g_l(matching_mass));
      
      %Replace RCM data with secondary experiment data (the primary 
      %experiment had issues with Agar precipitation)
      if strcmp(sample_media, 'RCM')
          mod_manifest.biomass(i) = ...
              secondary_biomass_data.wet_g_l(contains(secondary_biomass_data.media,sample_media));
      end
      
   end
end

end

