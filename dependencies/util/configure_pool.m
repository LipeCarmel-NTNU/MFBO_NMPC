function configure_pool(use_parallel, num_workers)
%CONFIGURE_POOL Open a process pool with num_workers, or ensure no pool exists.
%
%   configure_pool(true, n)  opens a pool with n workers. An existing pool
%   with a different worker count is deleted first.
%   configure_pool(false, ~) deletes any existing pool.

    p = gcp('nocreate');
    if use_parallel
        if isempty(p) || p.NumWorkers ~= num_workers
            if ~isempty(p)
                delete(p);
            end
            parpool('Processes', num_workers);
        end
    else
        if ~isempty(p)
            delete(p);
        end
    end
end
