function out = ra_tex_label(label)
%RA_TEX_LABEL Escape a controller label for the latex text interpreter.
%   ra_context sets the groot default text interpreter to latex, so a bare
%   "BO_1" in a title is read as a subscript command and errors. Only the
%   underscore needs escaping in these labels.
out = strrep(char(label), '_', '\_');
end
