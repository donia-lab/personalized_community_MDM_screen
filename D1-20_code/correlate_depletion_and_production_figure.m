%This script plots the depletion of drug against metabolite production and
%computes the corresponding correlation

%% Load prior metabolism analysis
clear;clc

metabolism_analysis = load('saved_analyses/drug_metabolism_analysis.mat');

include_all_metabolites = true;

%% Build table of production and depletion from original data

entropy_cutoff = 0.5;

%Import data
filtered_donor = metabolism_analysis.filtered_donor;
compound_list = metabolism_analysis.compound_list;

met_entropy = metabolism_analysis.met_entropy;
drug_entropy = metabolism_analysis.drug_entropy;
total_entropy = [met_entropy,drug_entropy];

%Make empty table
VariableNames = {'pearson','spearman','p_pearson','p_spearman','drug_vec','sum_met_vec',};
nan_matrix = nan(length(compound_list),length(VariableNames));
correlation_table = array2table(nan_matrix,'VariableNames',...
    VariableNames,'RowNames',compound_list);
num_donors = max(filtered_donor.donor);

template_vec = zeros(num_donors,1);
drug_vec_cell = {};
sum_met_vec_cell = {};

if include_all_metabolites
    met_labels = ...
        filtered_donor.Properties.VariableNames(contains(filtered_donor.Properties.VariableNames,'norm_met'));
else
    met_labels = {'norm_met_1'};
end

for i = 1:length(compound_list)
    
    compound = compound_list{i};
    compound_entropy = ...
        total_entropy(:,contains(total_entropy.Properties.VariableNames,compound));
    compound_donor = filtered_donor(contains(filtered_donor.compound,compound),:);
    drug_vec = template_vec;
    sum_met_vec = template_vec;
    
    for j = 1:num_donors
        donor_samples = compound_donor(compound_donor.donor == j,:);
        drug_vec(j) = mean(donor_samples.norm_drug);
        sum_met_vec(j) = mean(nansum(donor_samples{:,met_labels},2));
    end
    
    if max(compound_entropy.Variables) > entropy_cutoff
        [correlation_table{compound,'pearson'},correlation_table{compound,'p_pearson'}]=...
            corr(drug_vec,sum_met_vec);
        [correlation_table{compound,'spearman'},correlation_table{compound,'p_spearman'}]=...
            corr(drug_vec,sum_met_vec,'type','Spearman');
    end
    
    drug_vec_cell = [drug_vec_cell; drug_vec];
    sum_met_vec_cell = [sum_met_vec_cell; sum_met_vec];
    
end

correlation_table.drug_vec = drug_vec_cell;
correlation_table.sum_met_vec = sum_met_vec_cell;

filtered_table = correlation_table(~isnan(correlation_table.p_pearson),:);
filtered_table.corrected_p_pearson = ...
    mafdr(filtered_table.p_pearson,'BHFDR',true);
%% Plot production vs depletion

significance_cutoff = 0.01;
sig_vector = filtered_table.corrected_p_pearson < significance_cutoff;
bar_colors = [1 1 1; 0.7 0.7 0.7];

newfigure(4,2.5);

xlabels = filtered_table.Properties.RowNames;
hold on
for i = 1:length(sig_vector)
    h(i) = bar(i,filtered_table.pearson(i));
    h(i).FaceColor = bar_colors(sig_vector(i)+1,:);
end
ylabel({'Pearson correlation'})
xticks(1:length(sig_vector))
xticklabels(xlabels)
xtickangle(90)
box off
xlim([-0,length(xlabels)+1])
set(gca,'FontSize',9)
ylim([-1,1])
set(gca,'TickDir','out')

if include_all_metabolites
    print(gcf, '-dsvg', '-painters','figures/correlation_depletion_and_production_figure.svg');
else
    print(gcf, '-dsvg', '-painters','supp_figures/correlation_depletion_and_production_single_met_figure.svg');
end

%% Plot correlation coefficients
alpha = 0.5;
FontSize = 8;
color  = 'c';
for i = 1:length(filtered_table.Properties.RowNames)
    newfigure(1.5,1.5);
    compound = filtered_table.Properties.RowNames{i};
    drug_vec = correlation_table{compound,'drug_vec'}{1};
    sum_met_vec = correlation_table{compound,'sum_met_vec'}{1};
    
    y = sum_met_vec;
    x = [ones(size(drug_vec)),drug_vec];
    coeff = x\y;
    plot_x = sort(x);
    
    scatter(drug_vec,sum_met_vec,'o','MarkerEdgeColor','c','MarkerFaceColor',...
        'c','MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
    hold on
    plot(plot_x(:,2),plot_x*coeff,'k-')
    xlabel('Remaining drug','FontSize',8)
    ylabel('Metabolite','FontSize',8)
    ylims = get(gca,'Ylim');
    ylims(1) = 0;
    set(gca,'Ylim',ylims)
    set(gca,'FontSize',FontSize)
    title(compound,'FontSize',FontSize)
    set(gca,'TickDir','out')
    pause(3)
    if include_all_metabolites
        print(gcf, '-dsvg', '-painters',['figures/depletion_production_correlation/',compound,'.svg']);
    else
        print(gcf, '-dsvg', '-painters',['supp_figures/depletion_production_correlation_single_met/',compound,'_single_met.svg']);
    end
end


