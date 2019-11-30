%This script plots the results of the TPUP KO experiment and performs the
%relevant statistical testing. 

clear;clc
close all

measure = '250';

metab_peak = xlsread('TPUP_KO_data.xlsx',['metab_peak_',measure]);
parent_peak = xlsread('TPUP_KO_data.xlsx',['parent_peak_',measure]);
IS_peak = xlsread('TPUP_KO_data.xlsx',['IS_peak_',measure]);

%Data imported as TPUP UP TP WT going down

sum_metab = nansum(metab_peak,2);
sum_parent = nansum(parent_peak,2);
sum_IS = nansum(IS_peak,2); 

% final_norm = sum_metab./sum_IS;
% y_label = 'Normalized metabolite AUC';
% y_lim = 0.5;
% file_label = 'met_AUC';

final_norm = sum_metab./(sum_parent+sum_metab); 
y_label = 'Percent Conversion';
y_lim = 0.7;
file_label = 'conversion';

%Rearrange to WT TP UP TPUP horizontally
final_norm = [final_norm(10:12), final_norm(7:9), final_norm(4:6), final_norm(1:3)];

p(1) = 1;
[~,p(2)] = ttest2(final_norm(:,1),final_norm(:,2),'Vartype','unequal','tail','both');
[~,p(3)] = ttest2(final_norm(:,1),final_norm(:,3),'Vartype','unequal','tail','both');
[~,p(4)] = ttest2(final_norm(:,1),final_norm(:,4),'Vartype','unequal','tail','both');


mean_norm = nanmean(final_norm);
std_norm = nanstd(final_norm);

% Make plots and save
b = barwitherr(std_norm,mean_norm,'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0]);
set(gca,'box','off');
set(gca,'XTickLabels',{'WT','TP KO','UP KO','TPUP KO'})
ylabel(y_label);
xlim([0,5]);
%ylim([0,y_lim]);
set(gca,'FontSize',11)
set(gca,'TickDir','out');
set(gca,'FontName','SansSerif')

hold on 

space = 0.05;
t_space = 0.01;
bracket = 0.025;
tp = mean_norm(2); 
for i = 2:4
    y1 = tp+(i-1)*space;
    plot([1,i],[y1,y1],'k-')
    plot([1,1],[y1,y1-bracket],'k-')
    plot([i,i],[y1,y1-bracket],'k-')
    if p(i) >= 0.05
        text( (1+i)/2, y1 + 1.5*t_space, 'n.s.','HorizontalAlignment','center')
    elseif p(i) < 0.05 && p(i) >= 0.01
        text( (1+i)/2, y1 + t_space, '*','HorizontalAlignment','center')
    elseif p(i) < 0.01 && p(i) >= 0.001
        text( (1+i)/2, y1 + t_space, '**','HorizontalAlignment','center')
    elseif p(i) < 0.001
        text( (1+i)/2, y1 + t_space, '***','HorizontalAlignment','center')
    end
end

if strcmp(file_label, 'met_AUC')
    set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5]);
end

print(gcf, '-depsc2', ['TPUP_KO_', file_label,'_',measure,'.eps'])