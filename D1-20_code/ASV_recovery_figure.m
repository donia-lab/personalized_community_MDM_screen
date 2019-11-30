%Jaime G. Lopez 9/19. This figure computes the recovery of elements in the BG 
%cultures as a function of the abundance of the elements in the original feces.

%% Import data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%% Compute recoveries

BG_data = filtered_manifest(strcmp(filtered_manifest.media,'BG'),:);
feces_data = filtered_manifest(strcmp(filtered_manifest.media,'feces'),:);
num_samples = size(BG_data,1);
abundance_range = logspace(log10(1),log10(1e-4),101);
abundance_range(end+1) = 0;
fraction_entry_recovery = @(v1,v2,c) sum((v1 > 0)&(v2 > c))/sum(v2 > c);

mean_asv_recovery = nan(length(abundance_range),1);
mean_species_recovery = nan(length(abundance_range),1);
mean_genus_recovery = nan(length(abundance_range),1);
mean_family_recovery = nan(length(abundance_range),1);

for i = 1:length(abundance_range)
    asv_recovery = nan(num_samples,1);
    species_recovery = nan(num_samples,1);
    genus_recovery = nan(num_samples,1);
    family_recovery = nan(num_samples,1);
    for j = 1:num_samples
        donor = BG_data.donor(j);
        BG_asv = BG_data.rel_asv{j}.Variables; 
        feces_asv = feces_data(feces_data.donor == donor,:).rel_asv{1}.Variables;
        
        BG_species = BG_data.species{j}.Variables; 
        feces_species = feces_data(feces_data.donor == donor,:).species{1}.Variables;
        
        BG_genus = BG_data.genus{j}.Variables; 
        feces_genus = feces_data(feces_data.donor == donor,:).genus{1}.Variables;        
        
        BG_family = BG_data.family{j}.Variables; 
        feces_family = feces_data(feces_data.donor == donor,:).family{1}.Variables;             
      
        asv_recovery(j) = fraction_entry_recovery(BG_asv,feces_asv,abundance_range(i));
        species_recovery(j) = fraction_entry_recovery(BG_species,feces_species,abundance_range(i));
        genus_recovery(j) = fraction_entry_recovery(BG_genus,feces_genus,abundance_range(i));
        family_recovery(j) = fraction_entry_recovery(BG_family,feces_family,abundance_range(i));
    end
    mean_asv_recovery(i) = mean(asv_recovery);
    mean_species_recovery(i) = mean(species_recovery);
    mean_genus_recovery(i) = mean(genus_recovery);
    mean_family_recovery(i) = mean(family_recovery);
end


%% Plot recoveries

FontSize = 8;

newfigure(2,2);
hold on
LineWidth = 2;
plot(abundance_range(1:end-1),mean_family_recovery(1:end-1),'r-','LineWidth',LineWidth);
plot(abundance_range(1:end-1),mean_genus_recovery(1:end-1),'g-','LineWidth',LineWidth);
plot(abundance_range(1:end-1),mean_species_recovery(1:end-1),'b-','LineWidth',LineWidth);
plot(abundance_range(1:end-1),mean_asv_recovery(1:end-1),'k-','LineWidth',LineWidth);
xlabel({'Relative abundance', 'cutoff in feces'})
ylabel({'Mean fraction', 'recovered in BG'})
ylim([0,1])
set(gca,'XScale','log')
set(gca,'FontSize',FontSize)
leg = legend({'Family','Genus','Species','ASV'},'Location','southeast');
leg.ItemTokenSize = [10,18];
print(gcf, '-dsvg', '-painters','figures/ASV_recovery_figure.svg');

