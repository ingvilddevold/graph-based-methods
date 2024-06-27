function fig = plotHist(fig, H, varargin)
    % hist history from inner optimization
    opt = struct('legend', true);
    opt = merge_options(opt, varargin{:});

    if ~iscell(H)
        H = {H};
    end

    if ~ishandle(fig)
        fig = figure;
    else
        % Avoid stealing focus if figure already exists
        set(0, 'CurrentFigure', fig);
    end
    ax1 = subplot(2,1,1);
    title(ax1,'Objective');
    %set(ax1, 'XLim', xlim)
    set(ax1, 'YScale', 'log')
    hold(ax1,"on")
    grid(ax1,"on")

    ax2 = subplot(2,1,2);
    title(ax2, 'Gradient norm');
    %set(ax2, 'XLim', xlim)
    hold(ax2,"on")
    grid(ax2,"on")

    if opt.legend
        legend(ax1);
        legend(ax2);
    end

    popt = {'.-', 'LineWidth', 2};
    
    for i = 1:numel(H)
        hst = H{i};
        it = sum(hst.its, 'omitnan');

        xt = 0:(it-1);
        xlim = [-.2, xt(end)+.5];
        %ch = abs(hst.val(2:end)-hst.val(1:end-1));

        if isfield(hst, 'name') && ~isempty(hst.('name'))
            semilogy(ax1, xt, abs(hst.val(1:it)), popt{:}, 'DisplayName', hst.name);
            semilogy(ax2, xt, hst.pg(1:it), popt{:}, 'DisplayName', hst.name);
        else
            semilogy(ax1, xt, abs(hst.val(1:it)), popt{:});
            semilogy(ax2, xt, hst.pg(1:it), popt{:});
        end
    end

    % Add vertical lines representing graph refinements.
    for i = 1:(numel(H{1}.its)-1)
        x = cumsum(H{1}.its(1:i));
        xline(ax1, x-0.5, '--', 'Color', 'black', 'HandleVisibility', 'off')
        xline(ax2, x-0.5, '--', 'Color', 'black', 'HandleVisibility', 'off')
    end
    drawnow  
end
