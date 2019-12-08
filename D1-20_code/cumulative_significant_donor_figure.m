%This script makes a cumulative plot of total drugs vs number of
%significant donors

%% Import data from heatmap drug metabolism analysis

clear;clc

metabolism_analysis = load('saved_analyses/drug_metabolism_analysis.mat');
sig_met_table = metabolism_analysis.sig_met_table;
sig_depletion_table = metabolism_analysis.sig_depletion_table;

%% Loop through possible donor numbers

donor_num = 20;

cumulative_metabolites = zeros(donor_num+1,1);
cumulative_drugs = zeros(donor_num+1,1);

for i = 0:donor_num
    cumulative_metabolites(i+1) = sum(sum(sig_met_table.Variables) <= i);
    cumulative_drugs(i+1) = sum(sum(sig_depletion_table.Variables) <= i);
end

%% Plot cumulative curves

metabolite_y = reshape(repmat(cumulative_metabolites(:).',2,1),1,[]);
drug_y = reshape(repmat(cumulative_drugs(:).',2,1),1,[]);
x_val = 1:donor_num;
x_val = reshape(repmat(x_val(:).',2,1),1,[]);
x_val = [0,x_val,donor_num];

FontSize = 8;
LineWidth = 1.5;
newfigure(3,2);
drug_color = 'k';
met_color = 'r';
hold on

plot(x_val,drug_y/max(drug_y),'-','Color',drug_color,'LineWidth',LineWidth);
plot(x_val,metabolite_y/max(metabolite_y),'-','Color',met_color,'LineWidth',LineWidth);

ylabel('Fraction Cumulative Compounds','FontSize',FontSize)
xlabel('Number of significant donors','FontSize',FontSize)
set(gca,'FontSize',FontSize)

h = legend({'Drug depletion','Metabolite production'},'FontSize',FontSize,'Location','North');
h.Position(2) = 1.1*h.Position(2);
set(gca,'TickDir','out')
print(gcf, '-dsvg', '-painters','figures/cumulative_significant_donor_figure.svg');

