function filtered_manifest = filter_profinder_by_istd(manifest,compound_list)

%Identify internal standard
istd_properties = repmat([349.115,5.25],size(compound_list,1),1);
diff_from_istd = compound_list - istd_properties;
norm_from_istd = sqrt(diff_from_istd(:,1).^2 + diff_from_istd(:,2).^2);
istd_index = find(norm_from_istd == min(norm_from_istd));

%Normalize by internal standard
istd_cutoff = 5e5;
entity_indices = contains(manifest.Properties.VariableNames,'entity');
entities = manifest(:,entity_indices);
istd_areas = entities(:,['entity_',num2str(istd_index)]).Variables;
cutoff_failure_index = istd_areas < istd_cutoff;
istd_areas = repmat(istd_areas,[1,size(compound_list,1)]);
entities.Variables = entities.Variables./istd_areas;

manifest(:,entity_indices) = entities;

%Perform filtering
filtered_manifest = manifest(~cutoff_failure_index,:); 

end

