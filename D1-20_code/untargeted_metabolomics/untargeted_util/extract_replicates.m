function replicates = extract_replicates(manifest,donor,compound,entity)

entity_name = ['entity_',num2str(entity)];
compound_indices = strcmp(manifest.compound,compound);

if donor == 0
    replicates = manifest(compound_indices,entity_name).Variables;
else
    donor_indices = manifest.donor == donor;
    replicates = manifest(donor_indices & compound_indices,entity_name).Variables;
end

end

