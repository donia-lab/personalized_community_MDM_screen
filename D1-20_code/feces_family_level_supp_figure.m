%This script makes a family-level bar plot of all donor fecal samples

% Import data
clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);


%% Make stacked bar plots of samples by donor

other_cutoff = 0.05;

%bar_samples =  filtered_manifest(strcmp(filtered_manifest.media,'feces'),:);
bar_samples = sortrows(filtered_manifest,'donor','ascend');
%bar_samples =  bar_samples(contains(bar_samples.media,{'GAM','BG'}),:);
%bar_samples =  bar_samples(~(bar_samples.donor > 8 & contains(bar_samples.media,{'GAM'})),:);

family_matrix = [];

for i = 1:size(bar_samples,1)
    family_matrix = [family_matrix; transpose(bar_samples.family{i}.Variables)];
end


%Sort out families above cutoff and fully classified
orig_family_labels = bar_samples.family{i}.Properties.RowNames;
large_family_index = sum((family_matrix > other_cutoff)) > 0;
family_labels = orig_family_labels(large_family_index);
family_labels = split(family_labels,'f__');
family_labels = family_labels(:,2);
family_labels = family_labels(~strcmp(family_labels,''));

family_matrix = family_matrix(:,contains(orig_family_labels,family_labels));
family_matrix(:,end+1) = 1 - sum(family_matrix,2);

donor_list = unique(bar_samples.donor);

bracket_string = regexptranslate('wildcard','[*]');
contested = cellfun(@(x) regexp(x,bracket_string),family_labels,'UniformOutput',false);
for i = 1:length(family_labels)
    label = family_labels{i};
    if ~isempty(regexp(label,bracket_string,'ONCE'))
       family_labels{i} = label(2:end-1);
    end
end
family_labels = [family_labels;'Other'];
bar_colors = distinguishable_colors(length(family_labels));

%% Make actual bar plot

newfigure(6,4);

feces_index = strcmp(bar_samples.media,'feces');
feces_samples = bar_samples(feces_index,:);
feces_matrix = family_matrix(feces_index,:);

media_labels = strcat('D',num2str(feces_samples.donor));

h = bar(feces_matrix,'stacked');
box off
xticks(1:length(media_labels));
set(gca,'TickLabelInterpreter','none')
xticklabels(media_labels);
xtickangle(90)
xlim([0,length(media_labels)+1])
ylim([0,1])
yticks([0,1])
ylabel('Relative Abundance')

figure_size =  get(gcf, 'position');
h_legend = legend(family_labels,'Location','eastoutside');
set(h_legend, 'location', 'northeastoutside')


for j = 1:length(family_labels)
    h(j).FaceColor = 'flat';
    h(j).CData = bar_colors(j,:);
    h(j).EdgeColor = bar_colors(j,:);
end

set(gca,'FontSize',10)
h_legend.FontSize = 8;
pause(1)
print(gcf,'-dpng','supp_figures/feces_family_level_supp_figure.png','-r600');
