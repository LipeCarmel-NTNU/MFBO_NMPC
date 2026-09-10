function save_plot_outputs(figHandle, pngPath, fontSize, figWidthPx, figHeightPx)
%SAVE_PLOT_OUTPUTS Export standardized figure files (PNG 300 dpi + vector PDF).
arguments
    figHandle
    pngPath
    fontSize (1,1) double = 14
    figWidthPx (1,1) double = 900
    figHeightPx (1,1) double = 500
end
figure(figHandle);
set_fig_size(figWidthPx, figHeightPx);
set_font_size(fontSize);
exportgraphics(figHandle, pngPath, 'Resolution', 300);
[folderPath, fileStem] = fileparts(pngPath);
pdfPath = fullfile(folderPath, strcat(fileStem, '.pdf'));
save_figure(pdfPath, NaN, false);
end

%% ===================== GUARDED: REFINED FRONTIER =====================
