function [subject_table,modified_sample_table] = hmp_1_2_sample_to_subject(sample_table,coverage_cutoff)
%This function takes in the raw data for hmp-1-2 quantification and returns
%an aggregated table with unique subjects

%Get subject names from sample names and append to the sample table
subject_portions = cellfun(@(x) split(x,'-'),sample_table.sample,'UniformOutput',false);
subject_names = cellfun(@(x) x{9},subject_portions,'UniformOutput',false);
subject_SRS_ids = cellfun(@(x) x{10},subject_portions,'UniformOutput',false);
sample_table(:,'subject_id') = subject_names;
sample_table(:,'sample_id') = subject_SRS_ids;
sample_table.cohort = repmat({'HMP-1-2'},size(sample_table.gene));

subject_names = unique(subject_names);

%Make empty subject table
subject_labels = matlab.lang.makeValidName(subject_names);
subject_table = nan(length(subject_labels),4);
subject_table = array2table(subject_table,'RowNames',subject_labels,...
    'VariableNames',{'gene','coverage','RPKM','subject'});

for i = 1:length(subject_names)
    subject_name = subject_names{i};
    subject_label = subject_labels{i};
    
    %Get all samples related to subject
    subject_samples = sample_table(strcmp(sample_table.subject_id,subject_name),:);

    [coverage,RPKM] = aggregate_subject_data(subject_samples,coverage_cutoff);
    
    subject_table{subject_label,'coverage'} = coverage;
    subject_table{subject_label,'RPKM'} = RPKM;    
end

subject_table.subject = ...
                cellfun(@(x) ['hmp-1-2-stool-id-',strrep(x,'x','')],...
                subject_table.Properties.RowNames,'UniformOutput',false);

gene_cell = cell(size(subject_table,1),1);
gene_cell(:) = {sample_table.gene{1}};

subject_table.gene = gene_cell;
subject_table.Properties.RowNames = {};

%Rearrange sample table
modified_sample_table = sample_table(:,{'gene','coverage','RPKM','cohort','sample_id','subject_id'});

end