        plot(frBaseS(:,[2 2])', frBaseR(:,[1 3])', 'Color', .5*[1 1 1 .1]); hold on
        plot(frBaseS(:,[1 3])', frBaseR(:,[2 2])', 'Color', .5*[1 1 1 .1])

        plot(frBaseS(notSigBaseIx,2), frBaseR(notSigBaseIx,2), 'o', 'Color', [1 1 1], 'Linewidth', .25, 'MarkerSize', ms, 'MarkerFaceColor', [cmap(4,:)])
        plot(frBaseS(incBaseIx,2), frBaseR(incBaseIx,2), 'o', 'Color', [1 1 1], 'Linewidth', .25, 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
        plot(frBaseS(decBaseIx,2), frBaseR(decBaseIx,2), 'o', 'Color', [1 1 1], 'Linewidth', .25, 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))

        xlabel('Stationary Firing Rate')
        ylabel('Running Firing Rate')
         title('Base-driven firing rate')

        if strcmp(yscale, 'log')
            xlim([1 100])
            ylim([1 100])
            l = plot([1 100], [1 100]);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            l = plot([1 100], [1 100]*2);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            l = plot( 1:100, (1:100)/2);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            xlim([1 100])
            ylim([1 100])
        else
            xlim([0 100])
            ylim([0 100])
            l = plot([0 100], [0 100]);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            l = refline(1,20);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            l = refline(1,-20);
            set(l, 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
            xlim([0 100])
            ylim([0 100])
        end

        axis square
        plot.formatFig(gcf, [1.96 1.96], 'jnsci')
        set(gca, 'Yscale', yscale, 'Xscale', yscale)