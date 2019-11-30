%This figure looks at three simple model and computes the expected
%correlation between the known producer taxa and the resulting metabolite.

rt = 5;
single_producer = @(x) rt*x;
double_producer = @(x,y) rt*(x+y);
producer_consumer = @(x,z) ((x)./z).*(1-exp(-rt*z));
double_producer_consumer = @(x,y,z) ((x+y)./z).*(1-exp(-rt*z));

max_biomass = 1;
min_biomass = 0.3;

sample_num = 100;

rho1 = min_biomass + rand(sample_num,1)*(max_biomass-min_biomass);
rho2 = min_biomass + rand(sample_num,1)*(max_biomass-min_biomass);
rho3 = min_biomass + rand(sample_num,1)*(max_biomass-min_biomass);

data = [single_producer(rho1),double_producer(rho1,rho2),...
    producer_consumer(rho1,rho3),double_producer_consumer(rho1,rho2,rho3)];
my_alpha = 0.2;
my_color = 'r';
labels = {'(a)','(b)','(c)','(d)'};
newfigure(4,4);
for i = 1:size(data,2)
    subplot(2,2,i)
    scatter(rho1,data(:,i),'o','MarkerFaceColor',my_color,...
        'MarkerEdgeColor',my_color,'MarkerEdgeAlpha',my_alpha,'MarkerFaceAlpha',my_alpha)
    if sum(i == [3,4]) == 1
        xlabel('Producer biomass')
        xticks([0,1.5])
    else
        xticks([])
    end
    ylabel('Metabolite')
    xlim([0,1.5])
    ylim([0,round(1.2*max(data(:,i)))])
    yticks([0,round(1.2*max(data(:,i)))])
    set(gca,'FontSize',9)
    text(0.1,1,labels{i},'Units','normalized','FontWeight','bold')
    text(0.4,0.1,...
        ['\rho = ',num2str(round(corr(data(:,i),rho1,'type','Spearman'),2))],...
        'Units','normalized')
end

print(gcf, '-dpng','supp_figures/theoretical_taxa_vs_metabolism_supp_figure.png','-r600');


