clear;clc

%This script looks at the relationship between sequencing depth and ASV
%metrics in the PD sequencing data

%% Import and analyze

read_cutoff = 1e4;

%Read raw read count table
read_table = readtable('PD_read_counts.txt','ReadVariableNames',true,...
    'ReadRowNames',true);

%Read ASV table
asv_table = readtable('asv_table_split_read.txt','ReadRowNames',true);

%Select only main samples
asv_table = asv_table(:,contains(asv_table.Properties.VariableNames,'HD'));
asv_table = asv_table(:,~contains(asv_table.Properties.VariableNames,{'PBSr','FFC','PSBr'}));

%Preallocate vectors
depth = zeros(size(asv_table.Properties.VariableNames));
diversity = depth;

%Loop through samples and compute metrics
for i = 1:size(asv_table,2)
    sample_name = asv_table.Properties.VariableNames{i};
    sample = asv_table{:,sample_name};
    depth(i) = read_table{strrep(sample_name,'_','-'),'reads'};
    rel_sample = sample/sum(sample);
    diversity(i) = -nansum(rel_sample.*log2(rel_sample));
end

%Get media names from sample names
samples = asv_table.Properties.VariableNames;
media = cellfun(@(x) strsplit(x,'_'), samples,'UniformOutput',false);
media = cellfun(@(x) x{2}, media, 'UniformOutput',false);

%Turn data into table, filter out samples with very low reads
PD_table = table(media',diversity',depth','VariableNames',{'media','diversity','depth'});
PD_table = PD_table(PD_table.depth > 1e4,:);
PD_table = PD_table(~contains(PD_table.media,'Feces'),:);

%Compute correlation between diversity and depth
PD_depth_diversity_corr = corrcoef(PD_table.depth,PD_table.diversity);

%Perform one way anova of depth between media types
depth_anova = anova1(PD_table.depth,PD_table.media,'off');

%Fit LM of diversity including both media and depth
PD_lm = fitlm(PD_table,'diversity ~ 1 + depth + media');

PD_lm_anova = anova(PD_lm);
