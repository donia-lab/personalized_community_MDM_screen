%% Import data

clear;clc

manifest = wrapped_16S_import();

%% Identify samples that don't have enough reads to be skipped later

for i = 1:size(manifest,1)
    manifest.num_reads(i) =  sum(manifest.asv{i}.Variables);
    manifest.skip_sample(i) = manifest.num_reads(i) < 1e4;
end

%% Exclude replicates and blanks

rep_exclude = {'5.2','5.3','BG1','BG2','B.5.1','EMPTY','Blank'};

manifest = manifest(~contains(manifest.sample,rep_exclude),:);

%% Sort by media

sort_list = [];
media_order = {'feces','BB','BG','BHI','GAM','GMM','LB','liver','MRS','RCM','TB'};
for j = 1:max(manifest.donor)
    for i = 1:length(media_order)
        media_ind = strcmp(manifest.media,media_order{i});
        donor_ind = manifest.donor == j;
        sort_list = [sort_list; find(media_ind & donor_ind)];
    end
end

sorted_manifest = manifest(sort_list,:);

%% Generate asv level table and make clustergram for visualization

abundance_cutoff = 0.05;

asv_table = sorted_manifest.rel_asv{1,1};
for i = 2:size(sorted_manifest,1)
    asv_table = [asv_table, sorted_manifest.rel_asv{i,1}];
end
present_asv_index = sum(asv_table.Variables > abundance_cutoff,2) > 0;

filt_asv_table = asv_table(present_asv_index,:);

new_var = -log10(filt_asv_table.Variables);
new_var(isinf(new_var)) = 0;
f = clustergram(new_var,'Cluster', 'Column','Symmetric',false);
cluster_order = str2num(cell2mat(f.RowLabels));

%% Import taxonomy data

taxonomy = readtable('16S_data/taxonomy.csv','ReadRowNames',true,...
    'ReadVariableNames',true);


%% Pre-process data for heatmap
heatmap_data = filt_asv_table.Variables;
heatmap_data = heatmap_data(cluster_order,:);

%% Remove elements not identified at the order level

labels = filt_asv_table.Properties.RowNames(cluster_order);
asv_taxonomy = cell(length(labels),1);

for i = 1:length(labels)
    asv_taxonomy{i} = taxonomy{labels{i},'Taxon'}{1};
end
    
order_str = regexptranslate('wildcard','o__*; f');
order_labels  = cell(length(labels),1);
for i = 1:length(labels)
    [start_index,end_index] = regexp(asv_taxonomy{i},order_str);
    new_label = asv_taxonomy{i}(start_index+3:end_index-3);
    if isempty(start_index) || isempty(new_label)
        new_label = 'remove';
    end
    order_labels{i} =  new_label;
end

identified_orders = ~strcmp(order_labels,'remove');
heatmap_data = heatmap_data(identified_orders,:);
order_labels = order_labels(identified_orders);
labels = labels(identified_orders);

%% Add 'Other', log transform, and fix labels

order_labels{end+1} = 'Other';
labels{end+1} = 'Other';
heatmap_data(end+1,:) = 1 - nansum(heatmap_data,1);

heatmap_data(heatmap_data == 0) = NaN;
heatmap_data = log10(heatmap_data);

heatmap_data = flipud(heatmap_data);
order_labels = flipud(order_labels);
labels = flipud(labels);

media_order{strcmp(media_order,'liver')} = 'Liver';
media_order{strcmp(media_order,'feces')} = 'Feces';

labeled_heatmap_table = array2table(heatmap_data,'RowNames',labels,...
    'VariableNames',matlab.lang.makeValidName(sorted_manifest.sample));

%% Generate heatmap

color_max = log10(1);
color_min = log10(10^-5);
color_range = abs(color_max - color_min);
FontSize = 8;
pad_width = 3;
current_width = 0;
line_width = 1;
num_colors = 100;
colors = summer(num_colors);
value_to_index = @(x) round(num_colors*(x + color_range)/color_range);
%Make primary heatmap

newfigure(8,8/3);
hold on
for i = 1:max(manifest.donor)
    donor_data = heatmap_data(:,sorted_manifest.donor == i);
    skip_data = sorted_manifest.skip_sample(sorted_manifest.donor == i);
    for j = 1:size(donor_data,1)
        for k = 1:size(donor_data,2)
            if ~isnan(donor_data(j,k)) && ~skip_data(k)
                color = colors(value_to_index(donor_data(j,k)),:);
                rectangle('Position',[current_width + k, j+1, 1, 1],...
                    'FaceColor',color,'EdgeColor','none')
            end
        end
    end
    old_width = current_width;
    current_width = current_width + size(donor_data,2) + pad_width;
    label_position = (current_width -pad_width - 0.5*size(donor_data,2)+1);
    text(label_position,-10,['D',num2str(i)],'HorizontalAlignment','center',...
        'FontSize',FontSize)
    
    rectangle('Position',[old_width+1,1,size(donor_data,2),...
        size(donor_data,1)+2],'LineWidth',line_width)
    plot([label_position label_position],[1,-3],'k-')
end


%Add order level bar
order_list = unique(order_labels);
num_orders = length(order_list);
order_colors = flipud(distinguishable_colors(num_orders));

order_bar_width = 5;
order_bar_left_coord = -15;

for i = 1:length(order_labels)
    order_index = find(strcmp(order_list,order_labels{i}));
    rectangle('Position',[order_bar_left_coord, i+1, order_bar_width, 1],...
        'FaceColor',order_colors(order_index,:),'EdgeColor','none')
end
rectangle('Position',[order_bar_left_coord,1,...
    order_bar_width,length(order_labels)+2],'LineWidth',line_width);
text(order_bar_left_coord-5,1,'Taxonomic Order','Rotation',90,'FontSize',FontSize);

%Make legend for families
legend_left_coord = 0;
legend_top_coord = -25;
legend_spacing_y = 15;
legend_spacing_x = 65;
box_x = 5;
box_y = 5;
text_offset = 4;
box_width_x = 3;
box_width_y = 3;

tracking = 1;
for i = 1:box_x
    for j = 1:box_y
        if tracking <= length(order_list)
            element_x = legend_left_coord + (i-1)*legend_spacing_x;
            element_y =  legend_top_coord - legend_spacing_y*(j-1);
            rectangle('Position',[element_x,element_y - 0.5*box_width_y,...
                box_width_x,box_width_y],'EdgeColor','none','FaceColor',...
                order_colors(tracking,:))
            text(element_x + text_offset,element_y,order_list{tracking},...
                'HorizontalAlignment','left','VerticalAlignment','middle',...
                'FontSize',FontSize)
        end
        tracking = tracking + 1;
    end
end

%Add media label orders
text(1,size(donor_data,1)+15,...
    ['Sample order: ',strjoin(media_order,'|'),''],'FontSize',FontSize)

%Add colorbar
colorbar_left = 145;
colorbar_bottom = -50;
colorbar_height = 10;
for i = 1:num_colors
    rectangle('Position',[colorbar_left+i, colorbar_bottom, 1, colorbar_height],...
        'FaceColor',colors(i,:),'EdgeColor','none')
    
end
rectangle('Position',[colorbar_left+1, colorbar_bottom, size(colors,1), colorbar_height])
plot([colorbar_left+1 colorbar_left+1],[colorbar_bottom,colorbar_bottom-2],'k-')
text(colorbar_left+1,colorbar_bottom-12,'10^{-5}',...
    'HorizontalAlignment','center','FontSize',FontSize)
plot([colorbar_left+num_colors+1 colorbar_left+num_colors+1],[colorbar_bottom,colorbar_bottom-2],'k-')
text(colorbar_left+num_colors+1,colorbar_bottom-12,'10^{0}',...
    'HorizontalAlignment','center','FontSize',FontSize)

text(colorbar_left+1,colorbar_bottom + colorbar_height,'ASV Relative abundance',...
    'FontSize',FontSize,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom')

%Edit figure properties 
ylim([-90,150])
xlim([-20,current_width])
xticks([])
yticks([])
set(gca,'Visible','off')

print(gcf, '-dsvg', '-painters','figures/ASV_heatmap_figure.svg');
