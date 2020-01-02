%This script computes a version of the average ENDS figure with empirically
%estimated power

%% Make empirical estimate of the power curve

clear;clc

drug_analysis = load('saved_analyses/drug_metabolism_analysis.mat');
BG_table = drug_analysis.BG_table;
HK_table = drug_analysis.HK_table;
DMSO_table = drug_analysis.DMSO_table;

alpha = 0.01;

pass_list = [];
mean_list = [];
for i = 1:size(BG_table,1)
    for j = 1:drug_analysis.num_metabolites
        label = ['norm_met_',num2str(j),'_'];
        mean_label = [label,'mean'];
        p_label = [label,'p'];
        
        if ~isnan(BG_table{i,mean_label})
            pass_list = [pass_list; BG_table{i,p_label} < alpha; ...
                HK_table{i,p_label} < alpha;...
                DMSO_table{i,p_label} < alpha];
            mean_list = [mean_list; BG_table{i,mean_label};...
                HK_table{i,mean_label};...
                DMSO_table{i,mean_label}];
        end
    end
end

pass_list = pass_list(mean_list > 0);
mean_list = mean_list(mean_list > 0);

link = 'logit';
[b,~,stats] = glmfit(log(mean_list),pass_list,'binomial','link',link);
mean_range = logspace(-5,1,100);
pass_fit = glmval(b,log(mean_range),link);

figure
hold on
plot(mean_list,pass_list,'ko')
set(gca,'XScale','log')
plot(mean_range, pass_fit,'k-')

%% Compute average ENDS

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

% Define range of reaction rate values
r_lims = [-1,2];
r_list = logspace(r_lims(1),r_lims(2),40);

%Select only artificial media conditoins from data
media_list = unique(filtered_manifest.media);
media_list = media_list(~contains(media_list,{'feces','none'}));

%ENDS parameters and storage vector
alpha = 0.01;
n = 3;
mean_media_ENDS = zeros(length(media_list),length(r_list));
mean_best_media = cell(length(r_list),1);

empirical_ENDS = @(x,r) sum(glmval(b,log(r.*(x(x>0))),link));

for i = 1:length(r_list)
        
    %Compute local ENDS
    r = r_list(i);
    score_manifest = filtered_manifest;
    score_manifest.ENDS = zeros(size(score_manifest.biomass));
    for j = 1:size(score_manifest,1)
        absolute_asv = score_manifest.biomass(j).*score_manifest.rel_asv{j}.Variables;
        score_manifest.ENDS(j) = ...
            empirical_ENDS(absolute_asv,r);
    end
    % Find average ENDS per media
    for j = 1:length(media_list)
        media_cultures = score_manifest(strcmp(score_manifest.media,media_list{j}),:);
        mean_media_ENDS(j,i) = mean(media_cultures.ENDS);
    end
    
    mean_best_media{i} = media_list{mean_media_ENDS(:,i) == max(mean_media_ENDS(:,i))};
        
end

%% This part of the script loads the existing data and plots

FontSize = 8;
newfigure(3.5,2);
media_order = {'BB','BG','BHI','GMM','LB','liver','GAM','MRS','RCM','TB'};
media_labels = {'BB','BG','BHI','GMM','LB','Liver','mGAM','MRS','RCM','TB'};
colors = distinguishable_colors(length(media_list));
hold on

for i = 1:length(media_order)
    index = find(strcmp(media_list,media_order{i}));
    plot(r_list,mean_media_ENDS(index,:),'-','Color',colors(index,:),'LineWidth',2);
end
box off
yticks([0,20,40,60,80])
set(gca,'XScale','log')
l1 = legend(media_labels,'Location','eastoutside','FontSize',FontSize);
l1.ItemTokenSize = [10,18];
xlabel({'Reaction constant' ;'(normalized AUC per g/L biomass)'})
ylabel('Average ENDS')
xlim(10.^r_lims)
set(gca,'FontSize',FontSize)

print(gcf,'-dpng', '-painters','supp_figures/empirical_power_average_ENDS_supp_figure.png','-r600');

