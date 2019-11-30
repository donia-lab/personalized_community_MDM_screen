%% Import the pre-analyzed data from main heatmap script

clear;clc

metabolism_analysis = load('saved_analyses/drug_metabolism_analysis.mat');
sig_met_table = metabolism_analysis.sig_met_table;
filtered_donor = metabolism_analysis.filtered_donor;

donor_list = unique(filtered_donor.donor);

%% Make hydrocortisone metabolite figure

units = 'normalized AUC';

figure_x = 4*0.7;
figure_y = 1.7;
FontSize = 7;
FaceColors = {'none','r'};

newfigure(figure_x,figure_y);
hold on

compound = 'hydrocortisone';

hydro_samples = filtered_donor(strcmp(filtered_donor.compound,compound),:);
alpha = 0.5;

for donor = 1:length(donor_list)
    donor_label = ['donor_',num2str(donor)];
    donor_samples = hydro_samples(hydro_samples.donor == donor,:);
    donor_y = donor_samples.norm_met_1;
    donor_mean = mean(donor_y);
    donor_x = repmat(donor,size(donor_y));
    sig = sig_met_table{donor_label,[compound,'_met_1']};
    scatter(donor_x,donor_y,'o','MarkerEdgeColor','r','MarkerFaceColor',...
        FaceColors{sig+1},'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
    
end
xticks(1:max(donor_list))
set(gca,'TickLabelInterpreter','none')
xticklabels(strcat('donor_',string(1:20)))
xtickangle(90)
yticks([0,0.3,0.6])
ylim([0,0.6])
xlim([0.5,max(donor_list)+0.5]);
set(gca,'FontSize',FontSize);
ylabel({'Dihydrocortisone', units},'Interpreter','none')

print(gcf, '-dsvg', '-painters','figures/hydrocortisone_met_1_figure.svg');

%% Make sulfasalazine metabolite figure

newfigure(figure_x,figure_y);
hold on

compound = 'sulfasalazine';
FaceColors = {'none','b'};

sulfa_samples = filtered_donor(strcmp(filtered_donor.compound,compound),:);
alpha = 0.5;

for donor = 1:length(donor_list)
    donor_label = ['donor_',num2str(donor)];
    donor_samples = sulfa_samples(sulfa_samples.donor == donor,:);
    donor_y = donor_samples.norm_met_1;
    donor_mean = mean(donor_y);
    donor_x = repmat(donor,size(donor_y));
    sig = sig_met_table{donor_label,[compound,'_met_1']};
    scatter(donor_x,donor_y,'o','MarkerEdgeColor','b','MarkerFaceColor',...
        FaceColors{sig+1},'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
    
end
xticks(1:max(donor_list))
set(gca,'TickLabelInterpreter','none')
xticklabels(strcat('donor_',string(1:20)))
xtickangle(90)
yticks([0,0.5,1])
ylim([0,1])
xlim([0.5,max(donor_list)+0.5]);

set(gca,'FontSize',FontSize);

ylabel({'Sulfapyridine', units},'Interpreter','none')
%ylabel({'Hydrocortisone','met_1'})
print(gcf, '-dsvg', '-painters','figures/sulfasalazine_met_1_figure.svg');

%% Make capecitabine metabolite figure

newfigure(figure_x,figure_y);
hold on

compound = 'capecitabine';
FaceColors = {'none','g'};

cap_samples = filtered_donor(strcmp(filtered_donor.compound,compound),:);
alpha = 0.5;

for donor = 1:length(donor_list)
    donor_label = ['donor_',num2str(donor)];
    donor_samples = cap_samples(cap_samples.donor == donor,:);
    donor_y = donor_samples.norm_met_1;
    donor_mean = mean(donor_y);
    donor_x = repmat(donor,size(donor_y));
    sig = sig_met_table{donor_label,[compound,'_met_1']};
    scatter(donor_x,donor_y,'o','MarkerEdgeColor','g','MarkerFaceColor',...
        FaceColors{sig+1},'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)

end
xticks(1:max(donor_list))
set(gca,'TickLabelInterpreter','none')
xticklabels(strcat('donor_',string(1:20)))
xtickangle(90)
yticks([0,0.2,0.4])
ylim([0,0.4])
xlim([0.5,max(donor_list)+0.5]);

set(gca,'FontSize',FontSize);

ylabel({'Deglycocapecitabine', units},'Interpreter','none')

print(gcf, '-dsvg', '-painters','figures/capecitabine_met_1_figure.svg');


%% Make capecitabine metabolite figure


newfigure(figure_x,figure_y);
hold on

compound = 'digoxin';
FaceColors = {'none','m'};

cap_samples = filtered_donor(strcmp(filtered_donor.compound,compound),:);
alpha = 0.5;

for donor = 1:length(donor_list)
    donor_label = ['donor_',num2str(donor)];
    donor_samples = cap_samples(cap_samples.donor == donor,:);
    donor_y = donor_samples.norm_met_1;
    donor_mean = mean(donor_y);
    donor_x = repmat(donor,size(donor_y));
    sig = sig_met_table{donor_label,[compound,'_met_1']};
    scatter(donor_x,donor_y,'o','MarkerEdgeColor','m','MarkerFaceColor',...
        FaceColors{sig+1},'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)

end
xticks(1:max(donor_list))
set(gca,'TickLabelInterpreter','none')
xticklabels(strcat('donor_',string(1:20)))
xtickangle(90)
yticks([0,0.1,0.2])
ylim([0,0.2])
xlim([0.5,max(donor_list)+0.5]);

set(gca,'FontSize',FontSize);

ylabel({'Dihydrodigoxin', units},'Interpreter','none')

print(gcf, '-dsvg', '-painters','figures/digoxin_met_1_figure.svg');

