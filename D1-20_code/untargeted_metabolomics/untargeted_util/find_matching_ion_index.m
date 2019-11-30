function compound_index = find_matching_ion_index(target_ion,ion_list,mass_cutoff,rt_cutoff)
%This takes in a ion and ion list and find a match satisfying the cutoffs
%and also minimizing the norm of the vector of RT and mass differences.
%Returns empty if no matches passes cutoffs

num_ions = size(ion_list,1);

rep_rt = repmat([mass_cutoff,rt_cutoff],num_ions,1);
rep_ion = repmat(target_ion,num_ions,1);

ion_diff = abs(rep_ion - ion_list);
ion_norm  = sqrt(ion_diff(:,1).^2 + ion_diff(:,2).^2);

closest_ion = ion_norm == min(ion_norm);

passing_ions = sum(ion_diff < rep_rt,2) == 2;

compound_index = find(closest_ion & passing_ions);

end