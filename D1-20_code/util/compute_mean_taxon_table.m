function mean_taxon_table = ...
compute_mean_taxon_table(table,var,donor_order,truncate,prevalence_cutoff,abundance_cutoff)

%This function takes in a table with multiple BG replciates and returns a
%table of the mean composition of all taxon elements

%Make new taxon table
donor_names =  cellstr(strcat('donor_',string(donor_order)));
taxon_names = table{1,var}{1}.Properties.RowNames;
if truncate
    taxon_names = cellfun(@(x) x(round(length(x)/2):end),taxon_names,'UniformOutput',false);
end
taxon_names = matlab.lang.makeValidName(taxon_names);
taxon_names = matlab.lang.makeUniqueStrings(taxon_names);
mean_taxon_table = array2table(nan(length(donor_names),length(taxon_names)),'RowNames',...
    donor_names,'VariableNames',taxon_names);

%Add in taxa tables
for i = 1:length(donor_order)
    donor = donor_order(i);
    donor_samples = table(table.donor == donor,:);
    sum_composition = zeros(size(donor_samples{1,var}{1}.Variables));
    for j = 1:size(donor_samples,1)
        sum_composition = sum_composition +...
            donor_samples.biomass(j).*donor_samples{j,var}{1}.Variables;
    end
    mean_composition = transpose(sum_composition./size(donor_samples,1));
    mean_taxon_table(donor,:).Variables = mean_composition;
    
end

%Get only taxa present in more than samples than cutoff
prevalent_taxa = sum(mean_taxon_table.Variables>0,1) >= prevalence_cutoff;
mean_taxon_table = mean_taxon_table(:,prevalent_taxa);

abundant_taxa = sum(mean_taxon_table.Variables>abundance_cutoff,1) > 0;
mean_taxon_table = mean_taxon_table(:,abundant_taxa);

end

