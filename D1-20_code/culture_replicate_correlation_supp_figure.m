% This script looks at correlation in BG replicates at multiple levels

%% Import data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%% Loop through and compute correlations

media = 'GAM';

BG_manifest = filtered_manifest(strcmp(filtered_manifest.media,media),:);

compute_corr = @(x,y,z) corr(x(z,1).Variables,y(z,1).Variables,'Type','Pearson');

self_ASV = [];
self_shared = [];
nonself_ASV = [];
ASV_vec_1 = [];
ASV_vec_2 = [];
total_combos = nchoosek(1:size(BG_manifest,1),2);
for i = 1:size(total_combos,1)
    ind1 = total_combos(i,1);
    ind2 = total_combos(i,2);
    present_1 = BG_manifest.rel_asv{ind1}.Variables > 0;
    present_2 = BG_manifest.rel_asv{ind2}.Variables > 0;
    above_0_ASV = present_1 | present_2;
    
    ASV_corr = ...
        compute_corr(BG_manifest.rel_asv{ind1},BG_manifest.rel_asv{ind2},above_0_ASV);
    
    self = BG_manifest.donor(ind1) == BG_manifest.donor(ind2);
    if self
        self_ASV = [self_ASV; ASV_corr];
        self_shared = [self_shared; ...
            sum((present_1 & present_2))/sum(present_2); ...
            sum((present_1 & present_2))/sum(present_1)];
        ASV_vec_1 = [ASV_vec_1; BG_manifest.rel_asv{ind1}{above_0_ASV,:}];
        ASV_vec_2 = [ASV_vec_2; BG_manifest.rel_asv{ind2}{above_0_ASV,:}];
    else
        nonself_ASV = [nonself_ASV; ASV_corr];
    end
end


%% Look at correlations between BG and feces

feces_manifest = filtered_manifest(strcmp(filtered_manifest.media,'feces'),:);
feces_corr = [];
feces_shared = [];
for i = 1:size(feces_manifest,1)
    donor = feces_manifest.donor(i);
    partial_BG = BG_manifest(BG_manifest.donor == donor,:);
    present_1 = feces_manifest.rel_asv{i}.Variables > 0;
    for j = 1:size(partial_BG,1)
        present_2 = partial_BG.rel_asv{j}.Variables > 0;
        above_0_ASV = present_1 | present_2;
        ASV_corr = ...
            compute_corr(partial_BG.rel_asv{j},feces_manifest.rel_asv{i},above_0_ASV);
        feces_corr = [feces_corr; ASV_corr];
        feces_shared = [feces_shared; sum(present_1&present_2)];
    end
end


%% Plot overall relationship between one replicate and another

newfigure(3,3);
fig_alpha = 0.2;
scatter(ASV_vec_1,ASV_vec_2,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerEdgeAlpha',fig_alpha,'MarkerFaceAlpha',fig_alpha)
hold on
plot([1e-5,1],[1e-5,1],'k-','LineWidth',1.5)
set(gca,'YScale','log')
set(gca,'XScale','log')
text(1e-4,2*1e-1,['r = ',num2str(round(corr(ASV_vec_1,ASV_vec_2,'type','Pearson'),2))])
text(1e-4,2*10^(-1.5),['\rho = ',num2str(round(corr(ASV_vec_1,ASV_vec_2,'type','Spearman'),2))])
xlabel('Replicate 1 abundance')
ylabel('Replicate 2 abundance')
xticks([1e-5,1])
yticks([1e-5,1])
%title([media,' replicate correlation'])
set(gca,'FontSize',10)
print(gcf,'-dpng',['supp_figures/',media,'_replicate_correlation_figure.png'],'-r600');