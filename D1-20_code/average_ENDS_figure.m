%Jaime G. Lopez, 9/19. This script computes the average ENDS of all media
%across all samples and plots this as a function of the reaction rate

%% This portion of the script performs the analysis and saves it

% Import data
clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%Compute noise models for ENDS
donor_manifest_loc = 'QTOF_data/20x20_donor_manifest.xlsx';
donor_loc = 'QTOF_data/donor_data.csv';
istd = 'voriconazole';
plt = 0;
noise_models = estimate_noise_models(donor_manifest_loc,donor_loc,istd,plt);

% Define range of reaction rate values
r_lims = [-1,2];
r_list = logspace(r_lims(1),r_lims(2),40);

%Select only artificial media conditoins from data
media_list = unique(filtered_manifest.media);
media_list = media_list(~contains(media_list,{'feces','none'}));


%ENDS parameters and storage vector
alpha = 0.01;
n = 3;
only_ENDS = true;
read_cutoff = 0;
mean_media_ENDS = zeros(length(media_list),length(r_list));
mean_best_media = cell(length(r_list),1);

for i = 1:length(r_list)
        
    %Compute local ENDS
    r = r_list(i);
    score_manifest = ...
        compute_16S_scores(filtered_manifest,r,alpha,n,noise_models,read_cutoff,only_ENDS);
    
    % Find average ENDS per media
    for j = 1:length(media_list)
        media_cultures = score_manifest(strcmp(score_manifest.media,media_list{j}),:);
        mean_media_ENDS(j,i) = mean(media_cultures.ENDS);
    end
    
    mean_best_media{i} = media_list{mean_media_ENDS(:,i) == max(mean_media_ENDS(:,i))};
        
end

save('saved_analyses/average_ENDS_figure.mat','mean_media_ENDS','media_list','r_list','r_lims')

%% This part of the script loads the existing data and plots

load('saved_analyses/average_ENDS_figure.mat')
FontSize = 8;
newfigure(3.5,2);
media_order = {'BB','BG','BHI','GMM','LB','liver','GAM','MRS','RCM','TB'};
media_labels = {'BB','BG','BHI','GMM','LB','Liver','GAM','MRS','RCM','TB'};
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


print(gcf, '-dsvg', '-painters','figures/average_ENDS_figure.svg');

