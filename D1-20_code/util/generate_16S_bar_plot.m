function generate_16S_bar_plot(bar_samples,level,other_cutoff,donor,label_x,all_replicates,fontsize)

bar_samples = sortrows(bar_samples,'donor');

matrix = [];
for i = 1:size(bar_samples,1)
    matrix = [matrix; transpose(bar_samples{i,level}{1}.Variables)];
end

%Sort out families above cutoff and fully classified
orig_labels = bar_samples{i,level}{1}.Properties.RowNames;
large_index = sum((matrix > other_cutoff)) > 0;
labels = orig_labels(large_index);
labels = split(labels,[level(1),'__']);
labels = labels(:,2);
labels = labels(~strcmp(labels,''));

matrix = matrix(:,contains(orig_labels,labels));
matrix(:,end+1) = 1 - sum(matrix,2);

labels = [labels;'Other'];
bar_colors = distinguishable_colors(length(labels));

if all_replicates
    replicate_index = [];
    media_labels = {};
    bar_x = [];
    starting_x = 1;
    for i = 1:max(bar_samples.donor)
        index_bg = contains(bar_samples.media,{'BG'});
        index_feces = contains(bar_samples.media,{'feces'});
        index_d = bar_samples.donor == i;
        replicate_index = [replicate_index; find(index_feces & index_d);...
            find(index_bg & index_d)];
        media_labels = [media_labels ; strcat({'feces';'BG_1';'BG_2';'BG_3'},['_D',num2str(i)])];
        bar_x = [bar_x; (starting_x:(starting_x + 3))'];
        starting_x = starting_x + 5; 
    end
    bar_matrix = matrix(replicate_index,:);

else
    donor_index = bar_samples.donor == donor;
    bar_matrix = matrix(donor_index,:);
    bar_media = bar_samples.media(donor_index);
    media_labels = strcat(bar_media,['_D',num2str(donor)]);
    bar_x = 1:length(media_labels);
end

h = bar(bar_x,bar_matrix,'stacked');

box off

xticks(bar_x);

set(gca,'TickLabelInterpreter','none')
if label_x
    xticklabels(media_labels);
    xtickangle(90)
else
    xticklabels([])
end
xlim([0,max(bar_x)+1])
ylim([0,1])
yticks([0,1])
ylabel('Relative Abundance')

for j = 1:length(labels)
    h(j).FaceColor = 'flat';
    h(j).CData = bar_colors(j,:);
    h(j).EdgeColor = bar_colors(j,:);
end

h_leg = legend(labels,'Location','NorthEastOutside');
set(gca,'FontSize',fontsize)

h_leg.FontSize = 8;

end

