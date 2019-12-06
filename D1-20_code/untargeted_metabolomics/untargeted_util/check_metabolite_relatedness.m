function [link,met_present,drug_present] = ...
    check_metabolite_relatedness(drug_prop,met_prop,GNPS_data,cluster_info,mass_cutoff,rt_cutoff)
%This takes in a metabolite table, a parent molecule mass and RT, and a
%GNPS file and return whether the parent is present, the metabolite is
%present, and whether the two are related

cluster_ions = [cluster_info.x_ParentMass,cluster_info.x_RetTime/60];

drug_cluster_index = find_matching_ion_index(drug_prop,cluster_ions,mass_cutoff,rt_cutoff);
met_cluster_index = find_matching_ion_index(met_prop,cluster_ions,mass_cutoff,rt_cutoff);

drug_GNPS_index = [];
met_GNPS_index = [];

if ~isempty(drug_cluster_index)
    drug_GNPS_index = ...
        find(GNPS_data.clusterIndex == cluster_info.x_ClusterIdx(drug_cluster_index));
end
if ~isempty(met_cluster_index)
    met_GNPS_index = ...
        find(GNPS_data.clusterIndex == cluster_info.x_ClusterIdx(met_cluster_index));
end

if drug_GNPS_index == met_GNPS_index
    met_GNPS_index = [];
end
    
drug_present = ~isempty(drug_GNPS_index);
met_present = ~isempty(met_GNPS_index);

if drug_present && met_present
    drug_component = GNPS_data{drug_GNPS_index,'componentindex'};
    met_component = GNPS_data{met_GNPS_index,'componentindex'};
    drug_connected = drug_component ~= -1;
    
    link = (drug_component == met_component) & drug_connected; 
else
    link = false;
end

end


