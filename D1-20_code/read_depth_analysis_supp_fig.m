%This script looks at the relationship between sequencing depth and measured
%media properties in our screen

%% Import data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%% Add in the raw read numbers 

raw_read_file = '16S_data/raw_read_data.csv';
raw_reads = readtable(raw_read_file);
filtered_manifest.raw_reads = zeros(size(filtered_manifest.sample));
for i = 1:size(filtered_manifest,1)
   matching_index =  strcmp(raw_reads.sample,filtered_manifest.sample{i});
   filtered_manifest.raw_reads(i) = raw_reads.count(matching_index);
end

%% Compute ENDS and richness

%Compute noise models for ENDS
donor_manifest_loc = 'QTOF_data/20x20_donor_manifest.xlsx';
donor_loc = 'QTOF_data/donor_data.csv';
istd = 'voriconazole';
plt = 0;
noise_models = estimate_noise_models(donor_manifest_loc,donor_loc,istd,plt);

%ENDS parameters
alpha = 0.01;
n = 3;
only_ENDS = false;
read_cutoff = 0;
r = 0.5;

%Compute ENDS
score_manifest = ...
    compute_16S_scores(filtered_manifest,r,alpha,n,noise_models,read_cutoff,only_ENDS);


%% Make plot of ENDS vs depth

culture_manifest = score_manifest(~strcmp(score_manifest.media,'feces'),:);

newfigure(3,2);
color = 'g';
alpha = 0.2;
scatter(culture_manifest.raw_reads,culture_manifest.ENDS,'o',...
    'MarkerEdgeColor',color,'MarkerFaceColor',color,...
    'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
ylabel('ENDS')
xlabel('Reads sequenced')
text(3e3,70,['r = ',...
    num2str(round(corr(culture_manifest.raw_reads,culture_manifest.ENDS),2))]);
pause(1)
print(gcf, '-dpng','supp_figures/depth_vs_ENDS_supp_figure.png','-r600');

%% Make plot of richness vs depth

newfigure(3,2);
scatter(culture_manifest.raw_reads,culture_manifest.richness,'o',...
    'MarkerEdgeColor',color,'MarkerFaceColor',color,...
    'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
ylabel('ASV richness')
xlabel('Reads sequenced')
text(3e3,100,['r = ',...
    num2str(round(corr(culture_manifest.raw_reads,culture_manifest.richness),2))]);
pause(1)
print(gcf, '-dpng','supp_figures/depth_vs_richness_supp_figure.png','-r600');

%% Make jitter plot of read numbers
newfigure(3,2);
hold on
color = 'b';
alpha = 0.2;
jitter = @(x) x + 0.1*rand(size(x)) - 0.05; 
condition_list = unique(score_manifest.media);
for i = 1:length(condition_list)
    media_manifest = score_manifest(strcmp(score_manifest.media,condition_list{i}),:);
    ydata = media_manifest.raw_reads;
    scatter(jitter(repmat(i,size(ydata))),ydata,'o','MarkerEdgeColor',color,...
        'MarkerFaceColor',color,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
    plot([i-0.3,i+0.3],[mean(ydata),mean(ydata)],'k-','LineWidth',2)
end
xticks(1:length(condition_list))
xlim([0,length(condition_list)+1])
xticklabels(condition_list)
xtickangle(90)
ylim([0,6e4])
ylabel('Reads Sequenced')
set(gca,'FontSize',9)

print(gcf, '-dpng','supp_figures/read_count_jitter_plot_supp_figure.png','-r600');
