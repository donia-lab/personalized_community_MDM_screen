function [link,met_present,drug_present] = ...
    check_metabolite_relatedness(drug_prop,met_prop,GNPS_data,mass_cutoff,rt_cutoff)
%This takes in a metabolite table, a parent molecule mass and RT, and a
%GNPS file and return whether the parent is present, the metabolite is
%present, and whether the two are related

GNPS_ions = [GNPS_data.parentMass,GNPS_data.RTMean/60];

drug_index = find_matching_ion_index(drug_prop,GNPS_ions,mass_cutoff,rt_cutoff);
met_index = find_matching_ion_index(met_prop,GNPS_ions,mass_cutoff,rt_cutoff);

if drug_index == met_index
    met_index = [];
end
    
drug_present = ~isempty(drug_index);
met_present = ~isempty(met_index);

if drug_present && met_present
    drug_component = GNPS_data{drug_index,'componentindex'};
    met_component = GNPS_data{met_index,'componentindex'};
    drug_connected = drug_component ~= -1;
    
    link = (drug_component == met_component) & drug_connected; 
else
    link = false;
end

end


