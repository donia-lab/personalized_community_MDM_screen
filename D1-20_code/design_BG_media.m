% This script uses a prelimary set of communities (BB and GAM from donors
%1-9)in order to design a mixture such that diversity of the mixture of the
%two compositions is maximized. Note that we are using diversity here since
%BG was designed before ENDS was formalized.

%% Import data and select out only the original BB and GAM replicates

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

sample_codes = {'A.2','A.5.1'};

original_samples = filtered_manifest(contains(filtered_manifest.sample,sample_codes),:);

%% Find optimal mixtures of BB and GAM communities in all donors 

num_donors = 9; %BG designed based on sequencing of 9 initial donors 

combos = linspace(0,1,100);
shannon = @(P) nansum(-P.*log2(P));
opt_shannon_mix = zeros(1,num_donors);
shannon_array = zeros(length(combos),num_donors);
for donor = 1:num_donors
    
    donor_table = original_samples(original_samples.donor == donor,:); 
    
    BB_culture = donor_table(contains(donor_table.media,'BB'),:);
    BB_culture = BB_culture.rel_asv{1}.Variables;
    
    GAM_culture = donor_table(contains(donor_table.media,'GAM'),:);
    GAM_culture = GAM_culture.rel_asv{1}.Variables;

    for j = 1:length(combos)
        new_mix = combos(j)*BB_culture + (1-combos(j))*GAM_culture;
        shannon_array(j,donor) = shannon(new_mix); 
    end
    
    opt_shannon_mix(donor) = combos(shannon_array(:,donor) == max(shannon_array(:,donor)));
    
end


%% Make plot of diversity as a function of mixture

newfigure(3,3);
hold on
color = 'r';
for i = 1:size(shannon_array,2)
    plot(combos,shannon_array(:,i),color,'LineWidth',1.5);
end
plot([mean(opt_shannon_mix),mean(opt_shannon_mix)],[0,6],'k--','LineWidth',1.5)

xlim([0,1])
ylim([0,6])

xlabel('Fraction of BB')
ylabel('Shannon diversity')

print(gcf, '-dpng','supp_figures/BG_media_design_supp_figure.png','-r600');
