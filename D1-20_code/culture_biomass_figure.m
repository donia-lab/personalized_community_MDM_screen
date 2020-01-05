%Jaime G. Lopez, 9/19. This scripts makes boxplots of all of the culture
%biomasses by media type

%% Prep data for plotting

%Import data
biomass_file = '16S_data/primary_biomass_experiment.csv';
biomass_data = readtable(biomass_file);

biomass_list = [];
biomass_grp = [];

media_list = {'BB','BG','BHI','GMM','LB','liver','GAM','MRS','TB'};

%Aggregate into form usable by box plot function
for i = 1:length(media_list)
    culture_biomass = biomass_data.wet_g_l(strcmp(biomass_data.media,...
        media_list{i}));
    biomass_list = [biomass_list; culture_biomass];
    biomass_grp = [biomass_grp; zeros(length(culture_biomass),1) + i];
    
end

media_list{strcmp(media_list,'liver')} = 'Liver';
media_list{strcmp(media_list,'GAM')} = 'mGAM';

%% Plot
FontSize = 8;
newfigure(2,2);
old_position = get(gcf,'Position');
old_position(3) = 1.6*old_position(3);
set(gcf,'Position',old_position);
boxplot(biomass_list,biomass_grp,'Symbol','k+')
xticklabels(media_list)
xtickangle(90)
yticks([0,10,20,30,40])
ylim([0,40])
box off
ylabel('Biomass (g/L)')
box_line = 1;
set(gca,'FontSize',FontSize)

set(findobj(gcf,'tag','Median'), 'Color', 'k','LineWidth',box_line);
set(findobj(gcf,'tag','Box'), 'Color', 'k','LineWidth',box_line);
set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineWidth',box_line,'LineStyle','-');
set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineWidth',box_line,'LineStyle','-');
set(findobj(gcf,'tag','Upper Adjacent Value'), 'Color', 'k','LineWidth',box_line);
set(findobj(gcf,'tag','Lower Adjacent Value'), 'Color', 'k','LineWidth',box_line);

print(gcf, '-dsvg', '-painters','figures/culture_biomass_figure.svg');

