%This script generates family level plots for all donors and their
%corresponding cultures

%% Import the data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

filtered_manifest.media(contains(filtered_manifest.media,'GAM')) = {'mGAM'};

%% Make stacked bar plots of samples by donor

other_cutoff = 0.05;
tag = '';

bar_samples =  filtered_manifest(~strcmp(filtered_manifest.media,'none'),:);
bar_samples.media(strcmp(bar_samples.media,'liver')) = {'Liver'};
bar_samples = sortrows(bar_samples,'media');
feces_index = contains(bar_samples.media,'feces');
bar_samples = [bar_samples(feces_index,:); bar_samples(~feces_index,:)];

donor_list = unique(bar_samples.donor);

for i = 1:length(donor_list) 
    donor = donor_list(i);
    
    newfigure(6,8);
    [ha, pos] = tight_subplot(2,1,0.05,[0.13,0.05],[0.08,0.07]); 
    axes(ha(1));
    generate_16S_bar_plot(bar_samples,'order',other_cutoff,donor,false,false);
    axes(ha(2));
    generate_16S_bar_plot(bar_samples,'family',other_cutoff,donor,true,false);
    
    ax1_p = get(ha(1),'Position');
    ax2_p = get(ha(2),'Position');
    ax2_p(3:4) = ax1_p(3:4);
    set(ha(2),'Position',ax2_p);
    
    print(gcf,'-dpng',['supp_figures/family_and_order_level_plots/donor_',...
        num2str(donor),'_bar_plot',tag,'.png'],'-r600');
end


%% Make combined bar plot of all replicates

newfigure(16,6);

generate_16S_bar_plot(bar_samples,'family',other_cutoff,0,true,true);

print(gcf,'-dpng',['supp_figures/BG_replicate_figure.png'],'-r600');