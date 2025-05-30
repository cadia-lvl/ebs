function plotInPanel(ax, signal, gci, goi, tscl)

plot(ax, (0:length(signal)-1)*tscl, signal, 'k'); hold on;
yl = ylim;
y1 = 2*min(signal); 
y2 = 2*max(signal);
for ii=1:length(gci)-1
    patch([gci(ii);gci(ii);goi(ii);goi(ii)]*tscl,[y1; y2; y2; y1], 'k','LineStyle', 'none','FaceAlpha',0.3);
end
ylim(yl);
