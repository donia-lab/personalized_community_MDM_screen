function correlation_table = find_significant_correlations(coef_table,p_table,fdr)

p_list = p_table.Variables;

corrected_p_list = mafdr(p_list(:),'Showplot',false,'BHFDR','true');

corrected_p_table = p_table;
corrected_p_table.Variables = reshape(corrected_p_list,size(p_list));

item1 = {};
item2 = {};
sig_coef = [];
sig_p = [];


for i = 1:size(coef_table,1)
    for j = 1:size(coef_table,2)
        if corrected_p_table{i,j} < fdr
            item1 = [item1; coef_table.Properties.RowNames{i}];
            item2 = [item2; coef_table.Properties.VariableNames{j}];    
            sig_coef = [sig_coef; coef_table{i,j}];
            sig_p = [sig_p; corrected_p_table{i,j}];
        end
    end
end


correlation_table = table(item1,item2,sig_coef,sig_p);


end

