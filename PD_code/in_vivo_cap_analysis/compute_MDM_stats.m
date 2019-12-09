function [ob] = compute_MDM_stats(norm_data)

%This function computes statistics for the in vivo cap study

ob.av1 = nanmean(norm_data(1:6,:));
ob.av2 = nanmean(norm_data(7:end,:));

ob.dev1 = nanstd(norm_data(1:6,:));
ob.dev2 = nanstd(norm_data(7:end,:));


ob.n1 = [];
ob.n2 = [];
for i = 1:size(norm_data,2)
    ob.n1 = [ob.n1 sum(~isnan(norm_data(1:6,i)))];
    ob.n2 = [ob.n2 sum(~isnan(norm_data(7:end,i)))];
end

ob.sem1 = ob.dev1./sqrt(ob.n1);
ob.sem2 = ob.dev2./sqrt(ob.n2);

end

