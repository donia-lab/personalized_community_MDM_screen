%This script takes quantification results of genes on different microbiome
%cohorts and aggregates per cohort statistics

%% Set up tables and gene names

clear;clc

%Get all folders (each corresponds to gene)
all_files = dir('.');
gene_names = {'20a-HSDH-CLOSCI_00901_Clostridium-scindens-ATCC-35704',...
    '20b-HSDH-Bifido-HD-1-CL-100-DNA',...
    'deoA-BW25113_4382-E.coli-BW25113',...
    'udp-BW25113_3831-E.coli-BW25113'};
gene_labels = matlab.lang.makeValidName(gene_names);

%Specify the datasets that were quantified
dataset_names = {'diabetes-diseased','diabetes-healthy','fijicomp-Stool',...
    'hmp-1-2-Stool','hmp-stool','MetaHIT_Gut-danish-healthy',...
    'MetaHIT_Gut-spanish-CD','MetaHIT_Gut-spanish-healthy','MetaHIT_Gut-spanish-UC'};

gene_num = length(gene_labels);
dataset_num = length(dataset_names);

%Function to generate filenames
make_filename = @(gene,dataset) [gene,'/',gene,'-vs-',dataset,'-bowtie-quant.csv'];

%Make empty output tables
prevalence_table = nan(dataset_num,gene_num);
prevalence_table = array2table(prevalence_table,'RowNames',dataset_names,...
    'VariableNames',gene_labels);
abundance_table = prevalence_table;
summary_vars = ['N_total',cell(strcat('N_positive_',gene_names))];
summary_vars = matlab.lang.makeValidName(summary_vars);
summary_table = nan(dataset_num,length(summary_vars));
summary_table = array2table(summary_table,'VariableNames',summary_vars,....
    'RowNames',dataset_names);
full_data_table = cell(dataset_num,gene_num);
full_data_table = cell2table(full_data_table,'RowNames',dataset_names,...
    'VariableNames',gene_labels);

%% Loop through all data

coverage_cutoff = 0.5;
master_subject_table = [];
master_sample_table = [];
for i = 1:gene_num
    label = gene_labels{i};
    gene = gene_names{i};
    summary_var = summary_vars{i+1};
    for j = 1:dataset_num
        
        %Perform table import
        dataset = dataset_names{j};
        filename = make_filename(gene,dataset);
        sample_table = readtable(filename);
        sample_table.Properties.VariableNames = {'gene','coverage','RPKM','sample'};
        
        %Get table that have each row as a subject
        if strcmp(dataset,'hmp-stool')
            [subject_table,modified_sample_table] = hmp_1_1_sample_to_subject(sample_table,coverage_cutoff);
        elseif strcmp(dataset,'hmp-1-2-Stool')
            [subject_table,modified_sample_table] = hmp_1_2_sample_to_subject(sample_table,coverage_cutoff);
        else
            subject_table = sample_table;
            subject_table.Properties.VariableNames{4} = 'subject';
            modified_sample_table = sample_table;
            sample_strings = cellfun(@(x) split(x,'-'),sample_table.sample,'UniformOutput',false);
          
            %Make modified sample table
            if contains(dataset,'diabetes')
                modified_sample_table.cohort = ...
                    repmat({'Chinese'},size(modified_sample_table.gene));
                modified_sample_table.sample_id = ...
                    cellfun(@(x) x{4},sample_strings,'UniformOutput',false);
            elseif contains(dataset,'fiji')
                modified_sample_table.cohort = ...
                    repmat({'Fijicomp'},size(modified_sample_table.gene));
                modified_sample_table.sample_id = ...
                    cellfun(@(x) x{3},sample_strings,'UniformOutput',false);
            elseif contains(dataset,'danish')
                modified_sample_table.cohort = ...
                    repmat({'MetaHIT-Danish'},size(modified_sample_table.gene));
                modified_sample_table.sample_id = ...
                    cellfun(@(x) x{4},sample_strings,'UniformOutput',false);
            elseif contains(dataset,'spanish')
                modified_sample_table.cohort = ...
                    repmat({'MetaHIT-Spanish'},size(modified_sample_table.gene));
                modified_sample_table.sample_id = ...
                    cellfun(@(x) x{4},sample_strings,'UniformOutput',false);
            end
            
            modified_sample_table.subject_id = modified_sample_table.sample_id;
            modified_sample_table = ...
                modified_sample_table(:,{'gene','coverage','RPKM','cohort','sample_id','subject_id'});

        end
        
        coverage_index = subject_table.coverage > coverage_cutoff;
        summary_table{dataset,'N_total'} = size(subject_table,1);
        summary_table{dataset,summary_var} = sum(coverage_index);
        
        prevalence_table{dataset,label} = sum(coverage_index)/size(subject_table,1);
        abundance_table{dataset,label} = median(subject_table.RPKM(coverage_index));
        full_data_table{dataset,label} = {subject_table.RPKM(coverage_index)};
        
        master_subject_table = [master_subject_table;subject_table];
        master_sample_table = [master_sample_table;modified_sample_table];

    end
end
        
temp_var = abundance_table.Variables;
temp_var(isnan(temp_var)) = 0;
abundance_table.Variables = temp_var;

writetable(master_subject_table,'all_subject_table.csv','WriteRowNames',true,...
    'WriteVariableNames',true)

writetable(master_sample_table,'all_sample_table.csv','WriteRowNames',true,...
    'WriteVariableNames',true)

%% Aggregate the summary table into its final form

final_summary_table = summary_table('fijicomp-Stool',:);
final_summary_table.Properties.RowNames{...
    contains(final_summary_table.Properties.RowNames,'fijicomp-Stool')} = 'Fijicomp';
final_summary_table('HMP-1-1-Stool',:) = summary_table('hmp-stool',:);
final_summary_table('HMP-1-2-Stool',:) = summary_table('hmp-1-2-Stool',:);
final_summary_table{'Chinese',:} = sum(summary_table{{'diabetes-healthy','diabetes-diseased'},:});
final_summary_table{'MetaHIT-Spanish',:} = sum(...
    summary_table{{'MetaHIT_Gut-spanish-CD','MetaHIT_Gut-spanish-healthy','MetaHIT_Gut-spanish-UC'},:});
final_summary_table('MetaHIT-Danish',:) = summary_table('MetaHIT_Gut-danish-healthy',:);


writetable(final_summary_table,'summary_table.csv','WriteRowNames',true,...
    'WriteVariableNames',true)

%% Aggregate the data into its final form for heatmap plotting

countries = {'China','Fiji','USA1','USA2','Denmark','Spain'};
dataset_groups = {{'diabetes-diseased','diabetes-healthy'},'fijicomp-Stool',...
    'hmp-stool','hmp-1-2-Stool','MetaHIT_Gut-danish-healthy',...
    {'MetaHIT_Gut-spanish-CD','MetaHIT_Gut-spanish-healthy','MetaHIT_Gut-spanish-UC'}};

country_abundance = nan(length(countries),size(prevalence_table,2));
country_abundance = array2table(country_abundance,'RowNames',countries,...
    'VariableNames',prevalence_table.Properties.VariableNames);
country_prevalence = country_abundance; 

country_vec_table = cell(length(countries),size(prevalence_table,2));
country_vec_table = cell2table(country_vec_table,'RowNames',countries,...
    'VariableNames',prevalence_table.Properties.VariableNames);

for i = 1:size(country_abundance)
    country = countries{i};
    dataset_group = dataset_groups{i};
    country_prevalence_datasets = prevalence_table{dataset_group,:};
    country_abundance_datasets = abundance_table{dataset_group,:};
    country_full_data = full_data_table{dataset_group,:};
    
    country_temp_data = cell(1,gene_num);
    for j = 1:gene_num
        country_temp_data{j} = cat(1,country_full_data{:,j});
    end
    
    country_prevalence{country,:} = mean(country_prevalence_datasets,1);
    country_abundance{country,:} = mean(country_abundance_datasets,1);
    country_vec_table{country,:} = country_temp_data;
end


%% Make the capecitabine heatmap

cap_genes = {'deoA','udp'};
cap_index = contains(country_prevalence.Properties.VariableNames,cap_genes);
cap_prevalence = country_prevalence(:,cap_index);
cap_abundance = country_abundance(:,cap_index);
newfigure(3.5,2);
ha = tight_subplot(1, 2, 0.05, [0.1,0.05], [0.2,0.07]);
axes(ha(1));
imagesc(cap_prevalence.Variables);
colormap(gca,'parula')
h1 = colorbar('Location','northoutside');
h1.Limits = [0,1];
h1.Ticks = [0,1];
ylabel(h1,'Prevalence');
yticks(1:size(cap_abundance,1))
yticklabels(cap_abundance.Properties.RowNames);
xticklabels(cap_genes);
set(gca,'TickLength',[0,0])
hold on
size1 = size(cap_prevalence,1);
size2 = size(cap_prevalence,2);
for i = 1:size2+1
    plot([i-0.5,i-0.5],[0.5,size1+0.5],'k-')
    for j = 1:size1+1
        plot([0.5,size2+0.5],[j-0.5,j-0.5],'k-')
    end
end

axes(ha(2))
cmap_range_2 = [0,3];
imagesc(cap_abundance.Variables);
colormap(gca,'parula')
h2 = colorbar('Location','northoutside');
h2.Limits = cmap_range_2;
h2.Ticks = cmap_range_2;
yticks([])
xticklabels(cap_genes);
ylabel(h2,'Median RPKM');
pause(1)
set(gca,'TickLength',[0,0])
hold on
size1 = size(cap_prevalence,1);
size2 = size(cap_prevalence,2);
for i = 1:size2+1
    plot([i-0.5,i-0.5],[0.5,size1+0.5],'k-')
    for j = 1:size1+1
        plot([0.5,size2+0.5],[j-0.5,j-0.5],'k-')
    end
end
print(gcf, '-dsvg', '-painters','cap_genes_heatmap.svg');


%% Make the hydrocortisone heatmap

hydro_genes = {'20a','20b'};
hydro_index = contains(country_prevalence.Properties.VariableNames,hydro_genes);
hydro_prevalence = country_prevalence(:,hydro_index);
hydro_abundance = country_abundance(:,hydro_index);
newfigure(3.5,2);
ha = tight_subplot(1, 2, 0.05, [0.1,0.05], [0.2,0.07]);
axes(ha(1));
imagesc(hydro_prevalence.Variables);
colormap(gca,'parula')
h3 = colorbar('Location','northoutside');
h3.Limits = [0,1];
h3.Ticks = [0,1];
ylabel(h3,'Prevalence');
yticks(1:size(hydro_abundance,1))
yticklabels(hydro_abundance.Properties.RowNames);
xticklabels(hydro_genes);
set(gca,'TickLength',[0,0])
hold on
size1 = size(hydro_prevalence,1);
size2 = size(hydro_prevalence,2);
for i = 1:size2+1
    plot([i-0.5,i-0.5],[0.5,size1+0.5],'k-')
    for j = 1:size1+1
        plot([0.5,size2+0.5],[j-0.5,j-0.5],'k-')
    end
end

axes(ha(2))
imagesc(hydro_abundance.Variables);
colormap(gca,'parula')
h4 = colorbar('Location','northoutside');
h4.Limits = [0,2];
h4.Ticks = [0,2];
yticks([])
xticklabels(hydro_genes);
ylabel(h4,'Median RPKM');
pause(1)
hold on
size1 = size(hydro_prevalence,1);
size2 = size(hydro_prevalence,2);
for i = 1:size2+1
    plot([i-0.5,i-0.5],[0.5,size1+0.5],'k-')
    for j = 1:size1+1
        plot([0.5,size2+0.5],[j-0.5,j-0.5],'k-')
    end
end
print(gcf, '-dsvg', '-painters','hydro_genes_heatmap.svg');

%% Make jitter plots of cap abundance

cap_vec = country_vec_table(:,cap_index);

noise = 0.1;
jitter_x = @(len,x_pos) noise*rand(len,1) + repmat(x_pos,len,1) - 0.5*noise;

country_gap = 1;
gene_gap = 0.3;

newfigure(4,3);
hold on

colors = {'g','r'};
fig_alpha = 0.3;
current_x = 0;
xtick_set = zeros(size(cap_vec,1),1);
tracker = 1;
for i = 1:size(cap_vec,1)
    xtick_set(i) = current_x + 0.5*gene_gap;
    for j = 1:size(cap_vec,2)
        y_data = cap_vec{i,j}{1};
        plot_list(tracker) = ...
            scatter(jitter_x(length(y_data),current_x),y_data,'o',...
            'MarkerFaceColor',colors{j},'MarkerEdgeColor',colors{j},...
            'MarkerFaceAlpha',fig_alpha,'MarkerEdgeAlpha',fig_alpha);
        if j < size(cap_vec,2)
            current_x = current_x + gene_gap;
        end
        tracker = tracker + 1;
    end
    if i < size(cap_vec,1)
    current_x = current_x + country_gap;
    end
end

set(gca,'YScale','log')
xlim([0-0.5,current_x+0.5])
xticks(xtick_set)
xticklabels(countries)
xtickangle(90)
ylim([1e-2,1e3])
yticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
ylabel('RPKM')
legend(plot_list(1:2),cap_genes)
pause(1)
print(gcf, '-dsvg', '-painters','cap_genes_jitter_plot.svg');


%% Make jitter plots of hydro abundance

hydro_vec = country_vec_table(:,hydro_index);

newfigure(4,3);
hold on

current_x = 0;
xtick_set = zeros(size(hydro_vec,1),1);
tracker = 1;
colors = {'b','r'};
for i = 1:size(hydro_vec,1)
    xtick_set(i) = current_x + 0.5*gene_gap;
    for j = 1:size(hydro_vec,2)
        y_data = hydro_vec{i,j}{1};
        plot_list(tracker) = ...
            scatter(jitter_x(length(y_data),current_x),y_data,'o',...
            'MarkerFaceColor',colors{j},'MarkerEdgeColor',colors{j},...
            'MarkerFaceAlpha',fig_alpha,'MarkerEdgeAlpha',fig_alpha);
        if j < size(hydro_vec,2)
            current_x = current_x + gene_gap;
        end
        tracker = tracker + 1;
    end
    if i < size(hydro_vec,1)
    current_x = current_x + country_gap;
    end
end

set(gca,'YScale','log')
xlim([0-0.5,current_x+0.5])
xticks(xtick_set)
xticklabels(countries)
xtickangle(90)
ylim([1e-2,1e2])
yticks([1e-2,1e-1,1e0,1e1,1e2])
ylabel('RPKM')
legend(plot_list(1:2),hydro_genes,'Location','north')
pause(1)
print(gcf, '-dsvg', '-painters','hydro_genes_jitter_plot.svg');
