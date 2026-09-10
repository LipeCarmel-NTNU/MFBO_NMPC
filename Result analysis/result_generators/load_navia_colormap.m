function cmap = load_navia_colormap(n)
%LOAD_NAVIA_COLORMAP Load the sequential map used for fidelity-color plots.
matPath = which('navia.mat');
if strlength(string(matPath)) == 0
    error('Unable to locate navia.mat on MATLAB path.');
end
S = load(matPath, 'navia');
if ~isfield(S, 'navia')
    error('File %s does not contain variable ''navia''.', matPath);
end
cmap = double(S.navia);
if size(cmap, 2) ~= 3
    error('Variable ''navia'' in %s must have 3 columns (RGB).', matPath);
end
if size(cmap, 1) ~= n
    x  = linspace(0, 1, size(cmap, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(x, cmap, xq, 'linear');
end
cmap = min(max(cmap, 0), 1);
end
