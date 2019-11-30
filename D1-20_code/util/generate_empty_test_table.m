function [test_table,pval_types,ratio_types,mean_types,stat_info] = ...
    generate_empty_test_table(mol_types,donor_num,testable_compounds)

%This function generates and empty table for hypothesis testing

pval_types = strcat(mol_types,'_p');
ratio_types = strcat(mol_types,'_ratio');
mean_types = strcat(mol_types,'_mean');
stat_info = {'n_ctrl','n_exp'};
var_names = ['donor','drug',pval_types,ratio_types,mean_types,stat_info];
row_names = strcat(testable_compounds,['_',num2str(donor_num)]);
nan_matrix = nan(length(testable_compounds),length(var_names));

test_table = array2table(nan_matrix,'VariableNames',var_names,'RowNames',row_names);
test_table{:,'donor'} = donor_num;
test_table.drug = testable_compounds;


end

