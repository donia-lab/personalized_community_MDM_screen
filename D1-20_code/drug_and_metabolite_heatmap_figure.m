%% Import the data

clear;clc

donor_manifest_loc = 'QTOF_data/20x20_donor_manifest.xlsx';
donor_loc = 'QTOF_data/donor_data.csv';
BG_manifest_loc = 'QTOF_data/20x20_BG_manifest.xlsx';
BG_loc = 'QTOF_data/BG_data.csv';
HK_manifest_loc = 'QTOF_data/20x20_HK_manifest.xlsx';
HK_loc = 'QTOF_data/HK_data.csv';

pval_file = 'saved_analyses/drug_metabolism_analysis.mat';

istd = 'voriconazole';
num_metabolites = 12;

%Import QTOF quant data
donor_manifest = import_standard_QTOF_sample(donor_manifest_loc,donor_loc,istd,num_metabolites);
BG_manifest = import_standard_QTOF_sample(BG_manifest_loc,BG_loc,istd,num_metabolites);
HK_manifest = import_standard_QTOF_sample(HK_manifest_loc,HK_loc,istd,num_metabolites);
DMSO_manifest = import_DMSO_QTOF_sample(donor_manifest_loc,donor_loc,istd);

donor_list = unique(donor_manifest.donor);
compound_list = unique(donor_manifest.compound);
compound_list = compound_list(~cellfun('isempty',compound_list));

%% Add in dihydrodigoxin manual integration data

dihydrodigoxin_loc = 'QTOF_data/manual_digoxin_analysis.csv';
donor_manifest = add_dihydrodigoxin_data(donor_manifest,dihydrodigoxin_loc,0);
BG_manifest = add_dihydrodigoxin_data(BG_manifest,dihydrodigoxin_loc,0);
HK_manifest = add_dihydrodigoxin_data(HK_manifest,dihydrodigoxin_loc,0);
DMSO_manifest = add_dihydrodigoxin_data(DMSO_manifest,dihydrodigoxin_loc,1);


%% Remove outliers based on istd amount

filtered_donor = filter_outliers_from_istd(donor_manifest);
filtered_BG = filter_outliers_from_istd(BG_manifest);
filtered_HK = filter_outliers_from_istd(HK_manifest);
filtered_DMSO = filter_outliers_from_istd(DMSO_manifest);

%% Perform hypothesis testing donor by donor

BG_table = metabolism_hypothesis_testing(filtered_donor,filtered_BG,false,num_metabolites);
HK_table = metabolism_hypothesis_testing(filtered_donor,filtered_HK,false,num_metabolites);
DMSO_table = metabolism_hypothesis_testing(filtered_donor,filtered_DMSO,true,num_metabolites);


%% FDR correct using B-H

%Aggregate depletion pvals
depletion_pval = [BG_table.norm_drug_p; HK_table.norm_drug_p];

%Aggregate metabolite pvals
met_pval_str = regexptranslate('wildcard','norm_met_*_p');
met_pval_indices = regexp(BG_table.Properties.VariableNames,met_pval_str);
met_pval_indices = ~cellfun(@isempty,met_pval_indices);
met_pval_labels = BG_table.Properties.VariableNames(met_pval_indices);

metabolite_pval = [BG_table(:,met_pval_labels).Variables; ...
    HK_table(:,met_pval_labels).Variables;...
    DMSO_table(:,met_pval_labels).Variables];

%Correct pvalues for FDR
corrected_depletion_pval = mafdr(depletion_pval,'Showplot',false,'BHFDR','true');
corrected_metabolite_pval = mafdr(metabolite_pval(:),'Showplot',false,'BHFDR','true');
 
%Plot to see the correction
figure
subplot(1,2,1)
loglog(depletion_pval,corrected_depletion_pval,'ko')
xlabel('Original')
ylabel('Corrected')
title('Depletion p-values')
subplot(1,2,2)
loglog(metabolite_pval(:),corrected_metabolite_pval,'ko')
xlabel('Original')
ylabel('Corrected')
title('Metabolite p-values')
 
%Put the corrected pvalues into new tables
pval_length = size(BG_table,1);
BG_HK_labels = ['norm_drug_p',met_pval_labels];
corrected_BG_table = BG_table; 
corrected_BG_table(:,BG_HK_labels).Variables = nan(size(BG_table(:,BG_HK_labels)));
corrected_HK_table = HK_table; 
corrected_HK_table(:,BG_HK_labels).Variables = nan(size(HK_table(:,BG_HK_labels)));
corrected_DMSO_table = DMSO_table; 
corrected_DMSO_table(:,met_pval_labels).Variables = nan(size(DMSO_table(:,met_pval_labels)));

corrected_BG_table.norm_drug_p = corrected_depletion_pval(1:pval_length);
corrected_HK_table.norm_drug_p = corrected_depletion_pval(pval_length+1:end);

corrected_metabolite_pval = reshape(corrected_metabolite_pval,size(metabolite_pval)); 
corrected_BG_table(:,met_pval_labels).Variables = corrected_metabolite_pval(1:pval_length,:);
corrected_HK_table(:,met_pval_labels).Variables = corrected_metabolite_pval(pval_length+1:2*pval_length,:);
corrected_DMSO_table(:,met_pval_labels).Variables = corrected_metabolite_pval(2*pval_length+1:end,:);

%% Generate tables for heatmaps

max_donors = 20;

sig_pval_threshhold = 0.01;
sig_depletion_threshhold = 0.5;

nan_matrix = nan(max_donors,length(compound_list));

depletion_table = array2table(nan_matrix,'VariableNames',...
    compound_list,'RowNames',cellstr(strcat('donor_',string(fliplr(1:max_donors)))));
sig_depletion_table = depletion_table;

%Make metabolite labels for figure
metabolite_labels = cell(num_metabolites*length(compound_list),1);
counter = 1;
for i = 1:length(compound_list)
    for j = 1:num_metabolites
        metabolite_labels{counter+j-1} = [compound_list{i},'_met_',num2str(j)];
    end
    counter = counter + num_metabolites;
end


%Make empty metabolite heatmap table
met_nan_matrix = nan(max_donors,num_metabolites*length(compound_list));
mean_met_table = array2table(met_nan_matrix,'VariableNames',...
    metabolite_labels,'RowNames',cellstr(strcat('donor_',string(fliplr(1:max_donors)))));
sig_met_table = mean_met_table;

for i = 1:size(corrected_BG_table,1)
    drug = corrected_BG_table.drug{i};
    donor = ['donor_', num2str(corrected_BG_table.donor(i))];
    
    depletion_table{donor,drug} = corrected_BG_table.norm_drug_ratio(i);
    
    %Require drug to below depletion and pval threshold in all conditions
    sig_depletion_table{donor,drug} = ...
        (corrected_BG_table.norm_drug_ratio(i) < sig_depletion_threshhold) & ...
        (corrected_BG_table.norm_drug_p(i) <  sig_pval_threshhold) & ... 
        (corrected_HK_table.norm_drug_ratio(i) < sig_depletion_threshhold) & ...
        (corrected_HK_table.norm_drug_p(i) <  sig_pval_threshhold);
    for j = 1:num_metabolites
        metabolite = [drug,'_met_',num2str(j)];
        pval_mean_label = ['norm_met_',num2str(j),'_mean'];
        pval_label = ['norm_met_',num2str(j),'_p'];
        mean_met_table{donor,metabolite} = corrected_BG_table{i,pval_mean_label};
        
        %Require metabolite to be below pval threshold in all conditions
        sig_met_table{donor,metabolite} = ...
            (corrected_BG_table{i,pval_label} < sig_pval_threshhold) & ...
            (corrected_HK_table{i,pval_label} < sig_pval_threshhold) & ...
            (corrected_DMSO_table{i,pval_label} < sig_pval_threshhold);
    end
    
end

%Remove metabolites that aren't actually present
real_metabolite_index = ~isnan(mean_met_table(end,:).Variables);
mean_met_table = mean_met_table(:,real_metabolite_index);
sig_met_table = sig_met_table(:,real_metabolite_index);

%Normalize metabobolites to max
max_metabolites = repmat(max(mean_met_table.Variables),max_donors,1);
mean_met_table.Variables = mean_met_table.Variables./max_metabolites;

%Abbreviate mycophenolate mofetil name
depletion_table.Properties.VariableNames = ...
    strrep(depletion_table.Properties.VariableNames,...
    'mycophenolate_mofetil','myco_mofetil');
mean_met_table.Properties.VariableNames = ...
    strrep(mean_met_table.Properties.VariableNames,...
    'mycophenolate_mofetil','myco_mofetil');

%% Compute entropies for heatmaps
shannon = @(P) nansum(-P.*log2(P));

drug_entropy = depletion_table(1,:);
drug_entropy.Properties.RowNames = {'entropy'};
drug_entropy.Variables = zeros(1,size(drug_entropy,2));
for i = 1:size(sig_depletion_table,2)
    frac_significant = sum(sig_depletion_table(:,i).Variables)/size(depletion_table,1);
    drug_entropy{1,i} = shannon([frac_significant,1-frac_significant]);
end

met_entropy = mean_met_table(1,:);
met_entropy.Properties.RowNames = {'entropy'};
met_entropy.Variables = zeros(1,size(met_entropy,2));
for i = 1:size(sig_met_table,2)
    frac_significant = sum(sig_met_table(:,i).Variables)/size(mean_met_table,1);
    met_entropy{1,i} = shannon([frac_significant,1-frac_significant]);
end


save(pval_file);

%% Make drug depletion heatmap

FontSize = 8;

newfigure(4,4);
h = imagesc(depletion_table.Variables);
hold on
colormap parula
x_width = size(depletion_table,2);
y_width = size(depletion_table,1);
for i = 1:size(depletion_table.Variables,2)
    plot([i-0.5,i-0.5],[0.5,y_width+0.5],'k-')
    for j = 1:size(depletion_table.Variables,1)
        plot([0.5,x_width+0.5],[j-0.5,j-0.5],'k-')
        if sig_depletion_table{j,i} == 1
            plot(i,j,'k*','MarkerSize',3);
        end
    end
end
plot([i+1-0.5,i+1-0.5],[0.5,y_width+0.5],'k-')

%Make entropy axis
entropy_max = 3;
entropy_color = [0.7,0.7,0.7];
ax = gca;
box off
ax.Clipping = 'off';
entropy_zero = -0.5;
entropy_one = entropy_zero - entropy_max;
text_x = -0.25;
text(text_x,entropy_zero,'0','HorizontalAlignment','center','FontSize',FontSize);
plot([0.5,0.5],[entropy_zero,entropy_one],'k-')
text(text_x,entropy_one,'1','HorizontalAlignment','center','FontSize',FontSize);
text(text_x-2.5,mean([entropy_zero,entropy_one]),'Entropy',...
    'HorizontalAlignment','center','FontSize',FontSize);
for i = 1:size(sig_depletion_table,2)
    rectangle('Position',[i-0.5,entropy_zero-3*drug_entropy{1,i},1,3*drug_entropy{1,i}]...
        ,'EdgeColor','k','FaceColor',entropy_color)
end


set(gca,'TickLabelInterpreter','none')
set(gca,'XTick',1:size(depletion_table,2),'TickLength',[0,0])
set(gca,'XTickLabel',depletion_table.Properties.VariableNames)
set(gca,'YTick',1:size(depletion_table,1),'TickLength',[0,0])
set(gca,'YTickLabel',depletion_table.Properties.RowNames)
xtickangle(90)
set(gca,'TickLength',[0,0])

h2 = colorbar('eastoutside');
caxis([0,1])
ylabel(h2, 'Fraction drug remaining','FontSize',FontSize)

x1=get(gca,'position');
x=get(h2,'Position');
x(3)=0.04;
set(h2,'Position',x)
set(h2,'Ticks',[0,0.25,0.5,0.75,1])
set(h2,'TickLabels',{'0','0.25','0.5','0.75','1'},'FontSize',FontSize)
set(gca,'position',x1)

set(gca,'FontSize',FontSize)

pos=get(gca,'position'); 
pos(4)=0.9*pos(4);        
set(gca,'position',pos);  

print(gcf,'-dpng','figures/drug_heatmap_figure.png','-r600');
print(gcf, '-dsvg', '-painters','figures/drug_heatmap_figure_vector.svg');

%% Make metabolite heatmap

newfigure(1.2*4,4);
h1 = imagesc(mean_met_table.Variables);
hold on
colormap(parula)

set(gca,'TickLabelInterpreter','none')
set(gca,'XTick',1:size(mean_met_table,2))
set(gca,'XTickLabel',mean_met_table.Properties.VariableNames,'FontSize',FontSize)
set(gca,'YTick',1:size(mean_met_table,1))
set(gca,'YTickLabel',mean_met_table.Properties.RowNames,'FontSize',FontSize)
set(gca,'color',0*[1 1 1]);
xtickangle(90)
set(gca,'TickLength',[0,0])

x_width = size(mean_met_table,2);
y_width = size(mean_met_table,1);
for i = 1:size(mean_met_table.Variables,2)
    plot([i-0.5,i-0.5],[0.5,y_width+0.5],'k-')
    for j = 1:size(mean_met_table.Variables,1)
        plot([0.5,x_width+0.5],[j-0.5,j-0.5],'k-')
        if sig_met_table{j,i} == 1
            plot(i,j,'k*','MarkerSize',3);
        end
    end
end
plot([i+1-0.5,i+1-0.5],[0.5,y_width+0.5],'k-')

%Make entropy axis
ax = gca;
box off
ax.Clipping = 'off';
entropy_zero = -0.5;
entropy_one = entropy_zero - entropy_max;
text_x = -0.25;
text(text_x,entropy_zero,'0','HorizontalAlignment','center','FontSize',FontSize);
plot([0.5,0.5],[entropy_zero,entropy_one],'k-')
text(text_x,entropy_one,'1','HorizontalAlignment','center','FontSize',FontSize);
text(text_x-2.7,mean([entropy_zero,entropy_one]),'Entropy',...
    'HorizontalAlignment','center','FontSize',FontSize);
for i = 1:size(sig_met_table,2)
    rectangle('Position',[i-0.5,entropy_zero-3*met_entropy{1,i},1,3*met_entropy{1,i}]...
        ,'EdgeColor','k','FaceColor',entropy_color)
end


h12 = colorbar;
ylabel(h12,'Normalized AUC scaled to max','FontSize',FontSize)
x1=get(gca,'position');
x=get(h12,'Position');
x(3)=0.04;
set(h12,'Position',x)
set(h12,'Ticks',[0,0.25,0.5,0.75,1])
set(h12,'TickLabels',{'0','0.25','0.5','0.75','1'},'FontSize',FontSize)
set(gca,'position',x1)

set(gca,'FontSize',FontSize)

pos=get(gca,'position'); 
pos(4)=0.9*pos(4);        
set(gca,'position',pos);  


print(gcf,'-dpng','figures/metabolite_heatmap_figure.png','-r600');
print(gcf, '-dsvg', '-painters','figures/metabolite_heatmap_figure_vector.svg');

