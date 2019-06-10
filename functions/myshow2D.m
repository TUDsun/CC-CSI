function myshow2D(yinv, xinv, I, dBrange, fontsize, unitstr)

if isempty(dBrange)
    dBrange(1) = min(I(:));
    dBrange(2) = max(I(:));
    if dBrange(1) == dBrange(2)
        dBrange(2) = dBrange(2) + 1;
        dBrange(1) = dBrange(1) - 1;
    end
end

imagesc(yinv, xinv, I, dBrange)

xlabel(['$x$/' unitstr], 'interpreter', 'latex')
ylabel(['$y$/' unitstr], 'interpreter', 'latex')
axis equal tight
% axis xy
grid on
set(gca, 'layer', 'top')
% colormap(flipud(hot));
colormap jet
% caxis([0 2])
colorbar
% plottools('on')
% axis([min(yinv) +max(yinv) min(xinv) +max(xinv)])
set(gca, 'fontsize', fontsize, 'FontName', 'Times New Roman')
% ax = gca;
% ax.XMinorGrid = 'on';
% ax.YMinorGrid = 'on';
% ax.ZMinorGrid = 'on';