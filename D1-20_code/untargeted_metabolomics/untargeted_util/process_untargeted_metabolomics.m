function [putative_metabolite_table] = process_untargeted_metabolomics(donor_manifest_loc,...
    BG_manifest_loc,HK_manifest_loc,data_loc,donor_num,folder,load_existing_analysis)


if ~load_existing_analysis
    [donor_manifest,compound_list] = import_profinder(donor_manifest_loc,data_loc);
    BG_manifest = import_profinder(BG_manifest_loc,data_loc);
    HK_manifest = import_profinder(HK_manifest_loc,data_loc);
    
    % Add internal standard and filter outliers
    filtered_donor_manifest = filter_profinder_by_istd(donor_manifest,compound_list);
    filtered_BG_manifest = filter_profinder_by_istd(BG_manifest,compound_list);
    filtered_HK_manifest = filter_profinder_by_istd(HK_manifest,compound_list);
    % Hypothesis test all compounds in donor samples against BG, DMSO,and HK
    [BG_pval,HK_pval,DMSO_pval,BG_fc,HK_fc,DMSO_fc] = ...
        untargeted_hypothesis_testing(filtered_donor_manifest,...
        filtered_BG_manifest,filtered_HK_manifest,compound_list);
    save([folder,'D',num2str(donor_num),'_untargeted_metabolomics_pval.mat'])
else
    load([folder,'D',num2str(donor_num),'_untargeted_metabolomics_pval.mat']);
end

% Perform multiple hypothesis corrections
all_pval = [BG_pval(:,3:end);HK_pval(:,3:end);DMSO_pval(:,3:end)];
old_dim = size(all_pval);
all_pval = all_pval.Variables;

[~,all_qval] = mafdr(all_pval(:),'Showplot',false,'BHFDR','false');
all_qval = reshape(all_qval,old_dim);

BG_qval = BG_pval;
HK_qval = HK_pval;
DMSO_qval = DMSO_pval;

num_samples = size(BG_pval,1);
BG_qval(:,3:end).Variables = all_qval(1:num_samples,:);
HK_qval(:,3:end).Variables = all_qval(num_samples+1:2*num_samples,:);
DMSO_qval(:,3:end).Variables = all_qval(2*num_samples+1:3*num_samples,:);

% Make cutoffs based on log cutoff and qval

log2_cutoff = 1;
qval_cutoff = 0.01;

BG_sig = (BG_qval(:,3:end).Variables < qval_cutoff) & ...
    (BG_fc(:,3:end).Variables > log2_cutoff);

HK_sig = (HK_qval(:,3:end).Variables < qval_cutoff) & ...
    (HK_fc(:,3:end).Variables > log2_cutoff);

DMSO_sig = (DMSO_qval(:,3:end).Variables < qval_cutoff) & ...
    (DMSO_fc(:,3:end).Variables > log2_cutoff);

combined_sig = BG_sig & HK_sig & DMSO_sig;

%Exclude metabolites that are significant for multiple drugs
multidrug_metabolite = sum(combined_sig,1) > 1;
combined_sig(:,multidrug_metabolite) = 0; 

combined_sig_table = BG_pval;
combined_sig_table(:,3:end).Variables = combined_sig;

% Identify putative metabolites

drug_list = unique(BG_pval.compound);

putative_metabolite_table = array2table(NaN(length(drug_list),1),...
    'RowNames',drug_list,'VariableNames',{'metabolites'});

temp_cell = cell(length(drug_list),1);
for i = 1:length(drug_list)
    
    drug = drug_list{i};
    drug_index = strcmp(combined_sig_table.compound,drug);
    
    putative_index = sum(combined_sig_table(drug_index,3:end).Variables,1);
    putative_index = putative_index > 0;
    temp_cell{i} = compound_list(putative_index,:);
    
end

putative_metabolite_table.metabolites = temp_cell;

end
