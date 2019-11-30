%Jaime Lopez 9/19, this script generates the figure of the hydrocortisone
%gene discovery from the metagenomic library

%% Import all data
clear;clc

data_folder = 'QTOF_data/';
data_names = {'round_1_EF','round_1_GH','round_2','round_3','round_4','round_5_6'};
data_names = strcat(data_folder,data_names);
data_files = strcat(data_names,'.csv');
data_manifests = strcat(data_names,'_manifest.xlsx');
istd = 'voriconazole';

manifest = [];
for i = 1:length(data_files)
    partial_manifest = ...
        import_standard_QTOF_sample(data_manifests{i},data_files{i},istd);
    manifest = [manifest; sortrows(partial_manifest,'chosen','ascend')];
end

%% Plot the data on a jitter plot

newfigure(4,3);
hold on
cell_nums = [20000,2000,200,20,12,1];
lineages = {'E','F','G','H'};
round_manifest = manifest(contains(manifest.lineage,lineages),:);
clone_manifest = manifest(~contains(manifest.lineage,lineages),:);
colors = {'r','r','r','r'};
alpha = 0.3;
FontSize = 9;
for i = 1:size(round_manifest,1)
    met_level = round_manifest.norm_met_1(i);
    lineage = round_manifest.lineage{i};
    cell_num = round_manifest.cell_num(i);
    chosen = round_manifest.chosen(i);
    
    x_val = find(cell_nums == cell_num);
    MarkerEdgeColor = colors{strcmp(lineages,lineage)};
    
    jitter = (rand(1,1)-0.5)/4;
    
    if chosen
        scatter(x_val,met_level,'o','MarkerEdgeColor','k',...
            'MarkerFaceColor',MarkerEdgeColor,'MarkerEdgeAlpha',1,...
            'MarkerFaceAlpha',alpha);
    else
        scatter(x_val+jitter,met_level,'o','MarkerEdgeColor',MarkerEdgeColor,...
            'MarkerFaceColor',MarkerEdgeColor,'MarkerEdgeAlpha',alpha,...
            'MarkerFaceAlpha',alpha);
    end
end

offset = 0;

direct_met = clone_manifest.norm_met_1(strcmp(clone_manifest.lineage,'direct'));
scatter(x_val+offset+1,direct_met,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerEdgeAlpha',1,'MarkerFaceAlpha',alpha);
codon_met = clone_manifest.norm_met_1(strcmp(clone_manifest.lineage,'codon'));
scatter(x_val+offset+2,codon_met,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerEdgeAlpha',1,'MarkerFaceAlpha',alpha);
empty_met = clone_manifest.norm_met_1(strcmp(clone_manifest.lineage,'empty'));
scatter(x_val+offset+3,empty_met,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerEdgeAlpha',1,'MarkerFaceAlpha',alpha);

yticks([0,0.25,0.5,0.75])
ylim([0,0.75])
ticklabels = {'R1 (2x10^4)','R2 (2x10^3)','R3 (2x10^2)','R4 (2x10^1)',...
    'R5 (12)','R6 (1)','direct\_clone','codon\_opt','empty\_vec'};
xticklabels(ticklabels)
xticks(1:9)
xtickangle(90)
xlim([0.7,length(ticklabels)+0.3])

ylabel({'Hydrocortisone', 'metabolite signal'})
xlabel('Round type')

set(gca,'FontSize',FontSize)

print(gcf, '-dsvg', '-painters','figures/hydro_metagenomic_library_figure.svg');
