function ensure_dir(p)
%ENSURE_DIR Create a folder if it does not already exist.
    if ~isfolder(p)
        mkdir(p);
    end
end
