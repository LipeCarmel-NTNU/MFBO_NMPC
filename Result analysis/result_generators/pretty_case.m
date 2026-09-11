function label = pretty_case(name)
%PRETTY_CASE Folder name -> the arm label used in every figure and table.
%
%   results_cost_aware is the runtime-aware campaign. It was set up as a
%   multi-fidelity run, but the driver on the cluster carried no design-prefix
%   rows, so every training row sat at z = 1 and the runtime GP had no evidence
%   that cost varies with fidelity. All 38 of its evaluations ran at or beside
%   z = 1: what it actually demonstrates is cost-aware qLogNEHVI at full
%   fidelity, and the label says so rather than claiming multi-fidelity.
%
%   results_baseline is the runtime-unaware single-fidelity reference.
%
%   The folder names are what the registries and the results tree use, so the
%   mapping lives here rather than in a rename of anything the driver writes.
%   Anything else falls back to the old "Case n" / capitalised form, so a
%   campaign added later is still legible without editing this file.

    name = string(name);

    if contains(lower(name), "baseline")
        label = "SF";
        return
    end
    if contains(lower(name), "cost_aware")
        label = "Cost-aware";
        return
    end

    tok = regexp(name, 'case(\d+)\s*$', 'tokens', 'once');
    if isempty(tok)
        label = string(strrep(strrep(name, 'results_', ''), '_', ' '));
        if strlength(label) > 0
            label = upper(extractBefore(label, 2)) + extractAfter(label, 1);
        end
    else
        label = "Case " + string(tok{1});
    end
end
