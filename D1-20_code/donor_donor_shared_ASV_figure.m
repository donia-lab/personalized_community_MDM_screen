%This script makes a figure comparing the shared ASVs between donors and
%themselves and donors and other donors

%% Import data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%% Compute shared ASVs

feces_data = filtered_manifest(strcmp(filtered_manifest.media,'feces'),:);
num_donors = size(feces_data,1);

combos = nchoosek(1:num_donors,2);
donor_list = transpose(1:20);
combos = [combos;[donor_list,donor_list]];

[combined,grp,parametric_pval,exp_tstat] = compute_self_nonself_shared_ASVs(feces_data,feces_data,combos);

%% Compute permutation p value (samples are not independent!)

num_permutations = 3e3;

[permutation_pval,parametric_pval_array] = ...
    compute_permutation_pvalue(combined,combos,grp,num_permutations,exp_tstat);

%% Plot shared ASVs

FontSize = 8;

newfigure(1.5,2);
hold on
LineWidth = 1;

boxplot(combined,~grp,'Symbol','k+')
ylabel({'Shared ASVs'; '(Donor-Donor)'})
set(gca,'XTickLabel',{'Self','Non-self'})
ylim([0,1.05*max(combined)])
set(gca,'FontSize',FontSize)
box off

bar_offset = 10;
plot([1,2],[max(combined)+bar_offset, max(combined)+bar_offset],'k-');
plot([1,1],[max(combined)+bar_offset,max(combined)+bar_offset-3],'k-');
plot([2,2],[max(combined)+bar_offset,max(combined)+bar_offset-3],'k-');
ylim([0,max(combined)+bar_offset]);
text(1.5,max(combined)+bar_offset-2,'***','HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

set(findobj(gcf,'tag','Median'), 'Color', 'k','LineWidth',LineWidth);
set(findobj(gcf,'tag','Box'), 'Color', 'k','LineWidth',LineWidth);
set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineWidth',LineWidth,'LineStyle','-');
set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineWidth',LineWidth,'LineStyle','-');
set(findobj(gcf,'tag','Upper Adjacent Value'), 'Color', 'k','LineWidth',LineWidth);
set(findobj(gcf,'tag','Lower Adjacent Value'), 'Color', 'k','LineWidth',LineWidth);

ylim([0,200]);

xtickangle(90)

print(gcf, '-dsvg', '-painters','figures/donor_donor_shared_ASV_figure.svg');
