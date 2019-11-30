function [coefficient_table,p_table] = compute_table_correlations(tbl1,tbl2,corrtype)

%This function takes in tables of doubles and returns

var1 = tbl1.Properties.VariableNames; 
var2 = tbl2.Properties.VariableNames; 

nan_matrix = nan(length(var1),length(var2));
coefficient_table = array2table(nan_matrix, 'RowNames',var1,'VariableNames',var2);
p_table = coefficient_table; 

for i = 1:length(var1)
    for j = 1:length(var2)
        [coef,p] = corr(tbl1(:,i).Variables,tbl2(:,j).Variables,'Type',corrtype);
        coefficient_table{i,j} = coef;
        p_table{i,j} = p;
    end
end
end

