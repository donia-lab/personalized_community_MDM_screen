%This script generates the MN metadata for the file

%% Import necessary data and set up templates
clear;clc
hitting_sets = readtable('GNPS_analysis/hitting_donor_sets.csv');
abbrev_data = readtable('GNPS_analysis/drug_to_abbreviation.xlsx');

standard_template = 'BG-1X-P1-';

%% Generate the GROUP files

for i = 1:size(hitting_sets,1)
    drug = hitting_sets.drug_list{i};
    donors = hitting_sets(i,2:end).Variables;
    num_donors = sum(~isnan(donors));
    
    sample_loc = abbrev_data{strcmp(drug,abbrev_data.drug_list),'location'}{1};
    sample_abbrev = abbrev_data{strcmp(drug,abbrev_data.drug_list),'abbrev'}{1};

    standard_line = ['GROUP-1=',standard_template,sample_loc,'.mzXML'];
    group_new = {standard_line};
    
    second_line = 'GROUP-2=';
    for j = 1:num_donors
        new_sample = [sample_abbrev,'_D',num2str(donors(j)),'_',sample_loc,'.mzXML'];
        
        if j == 1
            second_line = [second_line,new_sample];
        else
            second_line = [second_line,',',new_sample];
        end

    end
    
    group_new = [group_new,second_line];
    
    filename = ['GNPS_analysis/GNPS_metadata/GROUP-',drug,'.txt'];
    
    filePh = fopen(filename,'w');
    fprintf(filePh,'%s\n',group_new{:});
    fclose(filePh);
    
end
