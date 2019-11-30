function [ENDS,bvec] = compute_ENDS(asv_table,r,alpha,n,noise_models,read_cutoff)
%This function computes the expected number of detectable strains for a given
%community. This assumes a power law error model and ignores singletons

%Remove singletons 
x = asv_table(asv_table > read_cutoff)./sum(asv_table);

%Power law std = a*mean^b
a = noise_models.a;
b = noise_models.b;

%Null distribution
null_mean = noise_models.null_mean;
null_std = noise_models.null_std;

est_std = @(xi) a*((r*xi).^b);

compute_detection_prob = @(xi) welch_power(null_mean,r*xi + null_mean,...
    null_std,sqrt(null_std^2 + est_std(xi)^2),n,alpha);

bvec = arrayfun(compute_detection_prob,x);

ENDS = sum(bvec); 
end

