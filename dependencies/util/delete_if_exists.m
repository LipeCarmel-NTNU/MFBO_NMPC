function delete_if_exists(p)
%DELETE_IF_EXISTS Delete a file if it exists. Missing files are ignored.
    if exist(p, 'file') == 2
        delete(p);
    end
end
