function format_tick(X_decimals, Y_decimals)
%FORMAT_TICK Format tick labels for linear and log axes.
%   - Linear axes: fixed decimal precision.
%   - Log axes   : labels as powers of 10 (or mantissa x 10^k).

    arguments
        X_decimals = []
        Y_decimals = []
    end

    ax = gca;
    holdState = ishold(ax);
    hold(ax, "on");

    if ~isempty(X_decimals)
        ax.XTickLabel = make_labels(ax.XTick, ax.XScale, X_decimals);
    end

    if ~isempty(Y_decimals)
        ax.YTickLabel = make_labels(ax.YTick, ax.YScale, Y_decimals);
    end

    if ~holdState
        hold(ax, "off");
    end
end


function labels = make_labels(ticks, scaleType, decimals)
    labels = strings(numel(ticks), 1);
    if strcmpi(scaleType, "log")
        for i = 1:numel(ticks)
            labels(i) = log_label(ticks(i), decimals);
        end
    else
        fmt = ['%.', num2str(decimals), 'f'];
        for i = 1:numel(ticks)
            labels(i) = string(sprintf(fmt, ticks(i)));
        end
    end
    labels = cellstr(labels);
end


function label = log_label(v, decimals)
    if ~(isfinite(v) && v > 0)
        label = "";
        return;
    end

    p = log10(v);
    pRound = round(p);

    if abs(p - pRound) < 1e-10
        label = "$10^{" + string(pRound) + "}$";
        return;
    end

    e = floor(p);
    m = v / 10^e;
    fmt = ['%.', num2str(decimals), 'f'];
    mStr = string(sprintf(fmt, m));
    label = "$" + mStr + "\times 10^{" + string(e) + "}$";
end
