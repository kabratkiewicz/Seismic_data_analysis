function H = Plot_Energy(S_energy, e_min, e_max, t_scale, f_scale, fontsize, max_size)

S = db(abs(S_energy)) - max(max(db(abs(S_energy))));
if (max_size == 1)
    H = figure('units','normalized','outerposition',[0 0 1 1]);
else
    H = figure;
end

imagesc(t_scale, f_scale, S);
set(gca,'YDir','normal')
xlabel('t [min]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [Hz]','FontSize',fontsize, 'interpreter','latex')
h=colorbar;title(h,'E[dB]','FontSize',fontsize, 'interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')
colormap(flipud(inferno(256)))
% colormap(viridis(256))
caxis([e_min e_max]) 
c = colorbar;
c.Label.String = 'E[dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
drawnow
end

