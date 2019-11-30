function [pearson,spearman] = plot_correlated_items(item1,item2,table1,table2)

x = table1{:,item1};
y = table2{:,item2};

figure
plot(x,y,'ko')
xlabel(item1,'Interpreter','none')
ylabel(item2,'Interpreter','none')

pearson = corr(x,y,'Type','Pearson');
spearman = corr(x,y,'Type','Spearman');

end

