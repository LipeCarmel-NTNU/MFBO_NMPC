function ensure_parent_dir(filepath)
%ENSURE_PARENT_DIR Create the folder that contains filepath, if needed.
    [parent, ~, ~] = fileparts(filepath);
    if parent ~= "" && ~isfolder(parent)
        mkdir(parent);
    end
end
