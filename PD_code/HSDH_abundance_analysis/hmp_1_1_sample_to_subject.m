function subject_table = hmp_1_1_sample_to_subject(sample_table,coverage_cutoff)
%This function takes in the raw data for hmp-1-1 quantification and returns
%an aggregated table with unique subjects

%Import the catalog and get only stool sample entries
catalog_file = 'hmp-1-1-project_catalog.csv';
catalog = readtable(catalog_file);
catalog = catalog(strcmp(catalog.HMPIsolationBodySite,'Gastrointestinal tract'),:);

%Make emtpy subject table
subject_names = unique(catalog.MRN_SubjectID_);
subject_labels = arrayfun(@num2str,subject_names,'UniformOutput',false);
subject_labels = matlab.lang.makeValidName(subject_labels);
subject_table = nan(length(subject_labels),4);
subject_table = array2table(subject_table,'RowNames',subject_labels,...
    'VariableNames',{'gene','coverage','RPKM','subject'});

%Append SRS IDs to the sample table for easier matching
sample_fragments = cellfun(@(x) split(x,'-'),sample_table.sample,'UniformOutput',false);
sample_table(:,'SRS_id') = cellfun(@(x) x{5},sample_fragments,'UniformOutput',false);

%Collect samples for each subject
for i = 1:length(subject_names)
    subject_name = subject_names(i);
    subject_label = subject_labels{i};
    
    %Get matching SRS IDs for the subject and use those to get samples
    subject_SRS_ids = ...
        catalog.SequenceReadArchiveID(catalog.MRN_SubjectID_ == subject_name);
    subject_samples = sample_table(contains(sample_table.SRS_id,subject_SRS_ids),:);
    
    %Return aggregated data for single subject
    [coverage,RPKM] = aggregate_subject_data(subject_samples,coverage_cutoff);
    
    subject_table{subject_label,'coverage'} = coverage;
    subject_table{subject_label,'RPKM'} = RPKM;    
end

subject_table = subject_table(~isnan(subject_table{:,1}),:);

subject_table.subject = ...
                cellfun(@(x) ['hmp-1-stool-id',strrep(x,'x','')],...
                subject_table.Properties.RowNames,'UniformOutput',false);

gene_cell = cell(size(subject_table,1),1);
gene_cell(:) = {sample_table.gene{1}};
subject_table.gene = gene_cell;
subject_table.Properties.RowNames = {};
end