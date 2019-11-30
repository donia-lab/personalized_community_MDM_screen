function fig = newfigure(width,height)

fig = figure;
set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf, 'PaperSize', [width height]);
set(gcf,'Position',[3,4,width,height]);
set(gcf,'PaperUnits','Inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperPosition',get(gcf,'Position'));
OuterPosition = get(gcf,'OuterPosition');
OuterPosition(4) = 0.9999*OuterPosition(4);
set(gcf,'OuterPosition',OuterPosition);

end