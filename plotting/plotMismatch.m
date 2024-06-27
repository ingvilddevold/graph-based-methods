function plotMismatch(H, varargin)
    % Mismatch plot for report.
    % Not including gradient
    
    if mod(nargin,2) == 0
        ax = varargin{1};
        varargin = varargin(2:end);
    else
        figure;
        ax = gca;
    end

    opt = struct('legend', true, ...
                 'verticalLines', false); % will use the first entry in H
    opt = merge_options(opt, varargin{:});

    if ~iscell(H)
        H = {H};
    end

    hold on
    box on

    %title(ax,'Mismatch vs. iterations');
    set(ax, 'YScale', 'log')

    if opt.legend
        legend(ax);
    end

    % Add vertical lines representing graph refinements.
    if opt.verticalLines
        for i = 1:(numel(H{1}.its)-1 )
            x = cumsum(H{1}.its(1:i));
            xline(ax, x-0.5, '--', 'Color', [0.7 0.7 0.7], ...
                'HandleVisibility', 'off', 'linewidth', 0.2)
        end
    end

    popt = {'.-', 'LineWidth', 2};

    maxVal = -Inf;
    minVal = Inf;
    maxIt = Inf;

    for i = 1:numel(H)
        hst = H{i};
        it = numel(hst.val);

        xt = 0:(it-1);
        maxVal = max(maxVal, max(H{i}.val));
        minVal = min(minVal, min(H{i}.val));
        maxIt = max(maxIt, it);

        if isfield(hst, 'name') && ~isempty(hst.('name'))
            semilogy(ax, xt, abs(hst.val(1:it)), popt{:}, 'DisplayName', hst.name);
        else
            semilogy(ax, xt, abs(hst.val(1:it)), popt{:});
        end
    end

    xlabel('Iterations')

    xlim = [-.2, maxIt];
    set(ax, 'XLim', xlim)
    set(ax, 'YLim', [minVal, maxVal])
end
