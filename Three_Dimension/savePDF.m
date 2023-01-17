% saving pdf figures
function savePDF(name)

h = gcf;
%set(gca,'FontSize',15)
box on
grid on

%set(gcf,'position',[1000,100,800,600])

set(h,'position',[100,100,900,600])
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
path = 'F:\MEGA\JP Collaboration\NEW codes HYPERVOLUME CONSERVE (june2021)\THREE DIM\Figure';
saveas(h,fullfile(path,[name,'.pdf'])) % saving fig in pdf
saveas(h,fullfile(path,[name,'.png'])) % saving fig in png

end
