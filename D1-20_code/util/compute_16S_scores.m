function mod_manifest = compute_16S_scores(manifest,r,alpha,n,noise_models,read_cutoff,only_ENDS)
%This function takes a manifest and computes various community scores

%Shannon entropy function
shannon = @(P) nansum(-P.*log2(P));

%JS divergence function
KLdiv = @(P,Q) nansum(P.*log(P./Q));
JSdiv = @(P,Q) 0.5*KLdiv(P,(Q+P)/2) + 0.5*KLdiv(Q,(Q+P)/2);

%Shared entry function
shared_entries = @(vec1,vec2) sum((vec1 > 0) & (vec2 > 0));

mod_manifest = manifest;
sample_names = manifest.sample;
num_sample = length(sample_names);

%Make empty new columns
manifest.ENDS = zeros(num_sample,1);
mod_manifest.richness_no_singleton = zeros(num_sample,1);

if ~only_ENDS
    mod_manifest.shannon = zeros(num_sample,1);
    mod_manifest.richness = zeros(num_sample,1);
    mod_manifest.shared_asv = zeros(num_sample,1);
    mod_manifest.JS = zeros(num_sample,1);
    mod_manifest.shared_asv_above_1_pct = zeros(num_sample,1);
    mod_manifest.b_vector = cell(num_sample,1);
end

%Loop through and compute scores per sample
for i = 1:num_sample
    rel_table = manifest{i,'rel_asv'}{1};
    rel_table = rel_table.Variables;
    rel_family = manifest{i,'family'}{1};
    rel_family = rel_family.Variables;
    count_table = manifest{i,'asv'}{1};
    count_table = count_table.Variables;
    biomass = manifest{i,'biomass'};
    
    [mod_manifest.ENDS(i), mod_manifest.b_vector{i}] = ...
        compute_ENDS(count_table,biomass.*r,alpha,n,noise_models,read_cutoff);
    mod_manifest.richness_no_singleton(i) = sum(count_table > 1);
    
    if ~only_ENDS
        mod_manifest.shannon(i) = shannon(rel_table);
        mod_manifest.richness(i) = sum(count_table > 0);
        if manifest.donor(i) > 0
            feces_loc = contains(manifest.media,'feces') ...
                & (manifest.donor == manifest.donor(i));
            rel_feces = manifest.rel_asv{feces_loc};
            rel_feces = rel_feces.Variables;
            rel_feces_family = manifest.family{feces_loc};
            rel_feces_family = rel_feces_family.Variables;      
            mod_manifest.JS(i) = JSdiv(rel_family,rel_feces_family);
            mod_manifest.shared_asv(i) = shared_entries(rel_table,rel_feces);
            mod_manifest.shared_asv_pct(i) = mod_manifest.shared_asv(i)/sum(rel_feces > 0);
            mod_manifest.shared_asv_above_1_pct(i) = ...
                shared_entries(rel_table,rel_feces.*(rel_feces>0.01))/sum(rel_feces > 0.01);
        end
    end
end

end

