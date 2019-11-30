function [coverage,RPKM] = aggregate_subject_data(subject_samples,coverage_cutoff)
%This function takes samples from a single donor and returns an aggregated
%coverage and RPKM. If there is only one same present, that sample's data
%is returned. If there is more than one sample and none of them meet the 
%coverage threshold, then the average of samples is shown. If there
%is/are sample(s) above the coverage threshold, their average is
%reported

coverage_vec = subject_samples.coverage;
RPKM_vec = subject_samples.RPKM;

coverage_index = coverage_vec > coverage_cutoff;

if size(subject_samples,1) == 1 || sum(coverage_index) == 0
    coverage = mean(coverage_vec);
    RPKM = mean(RPKM_vec);
elseif sum(coverage_index) > 0
    coverage = mean(coverage_vec(coverage_index));
    RPKM = mean(RPKM_vec(coverage_index));
end


end

