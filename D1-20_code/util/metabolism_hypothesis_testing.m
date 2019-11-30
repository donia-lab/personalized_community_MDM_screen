function testing_table = metabolism_hypothesis_testing(donor_manifest,ctrl_manifest,DMSO,num_met)

%This function performs hypothesis testing for MDM drug metabolism

testing_table = [];

donor_list = unique(donor_manifest.donor);

mol_types = cellstr(strcat('norm_met_',string(1:num_met)));
tail_types = cell(size(mol_types));
tail_types(:) = {'right'};

%Define types of molecules and the type of ttest
if ~DMSO
    mol_types = ['norm_drug',mol_types];
    tail_types = ['left',tail_types];
end


for i = 1:length(donor_list)
    
    donor_num = donor_list(i);
    
    %Get data for one donor
    single_donor = donor_manifest(donor_manifest.donor == donor_num,:);
    donor_set = single_donor.set(1);
    if DMSO
        single_ctrl = ctrl_manifest(ctrl_manifest.donor == donor_num,:);
    else
        single_ctrl = ctrl_manifest(ctrl_manifest.set == donor_set,:);
    end
    
    
    %List of compounds that are in both experimental and control
    testable_compounds = unique(single_donor.compound);
    
    [partial_table,pval_types,ratio_types,mean_types,stat_info] = ...
        generate_empty_test_table(mol_types,donor_num,testable_compounds);
    
    row_names = partial_table.Properties.RowNames;
    
    %Fill table with p values and relevant values
    for j = 1:length(row_names)
        for k = 1:length(mol_types)
            compound_index = row_names{j};
            compound = split(compound_index,'_');
            compound = strjoin(compound(1:end-1),'_');
            donor_sample = ...
                single_donor(strcmp(single_donor.compound,compound),mol_types{k}).Variables;
            
            if DMSO
                metab_num = split(mol_types{k},'_met_');
                met_name = ['norm_',compound,'_met_',metab_num{end}];
                ctrl_sample = ...
                    single_ctrl(:,strcmp(single_ctrl.Properties.VariableNames,met_name)).Variables;
            else
                ctrl_sample = ...
                    single_ctrl(strcmp(single_ctrl.compound,compound),mol_types{k}).Variables;
            end
            
            if ~isempty(ctrl_sample)
                [~,partial_table{compound_index,pval_types{k}}] = ...
                    ttest2(donor_sample,ctrl_sample,'Vartype','unequal','tail',tail_types{k});
                partial_table{compound_index,ratio_types{k}} = mean(donor_sample)/mean(ctrl_sample);
                partial_table{compound_index,mean_types{k}} = mean(donor_sample);
                partial_table{compound_index,stat_info{1}} = size(ctrl_sample,1);
                partial_table{compound_index,stat_info{2}} = size(donor_sample,1);
            end
        end
    end
    
    testing_table = [testing_table; partial_table];
    
end



