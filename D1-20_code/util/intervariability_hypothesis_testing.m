function intervariability_table = intervariability_hypothesis_testing(donor_manifest)

%This script tests for intervariability in metabolite and drug depletion
%between different donors

compound_list = unique(donor_manifest.compound);
donor_list = unique(donor_manifest.donor);

mol_types = {'norm_drug','norm_met_1','norm_met_2'};
table_vars = strcat(mol_types,'_p');
table_rows = compound_list;
nan_matrix = nan(length(table_rows),length(table_vars));
intervariability_table = array2table(nan_matrix,'VariableNames',...
    table_vars,'RowNames',table_rows);

for i = 1:length(compound_list)
    compound = compound_list{i};
    compound_data = donor_manifest(strcmp(donor_manifest.compound,compound),:);
    for type = 1:length(mol_types)  
        sample_list = [];
        label_list = [];
        for donor = 1:max(donor_list)
            donor_data = compound_data(compound_data.donor == donor,:);
            type_data = donor_data{:,mol_types{type}};
            sample_list = [sample_list; type_data];
            label_list = [label_list; zeros(size(type_data))+donor];
        end
        intervariability_table{compound,table_vars{type}} = ...
            anova1(sample_list,label_list,'off');
    end
end

end

