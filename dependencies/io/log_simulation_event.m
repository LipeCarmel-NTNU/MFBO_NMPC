function log_simulation_event(log_path, message)
%LOG_SIMULATION_EVENT Append a line to a log file without ever aborting the run.
%
%   log_simulation_event(log_path, message) appends message as one line.
%   log_simulation_event(log_path) retries any lines still buffered from
%   earlier failed writes, which is useful once at the end of a run.
%
%   The log is a diagnostic side channel, so a failure to write must not
%   interrupt a simulation that may already have run for hours. On Windows an
%   editor holding the file open makes fopen fail, so this function retries a
%   few times, then keeps the line in memory and writes it on a later call.
%   Nothing here throws: every failure path ends in a warning at most, and
%   the warning is issued once per session so a locked file cannot flood the
%   command window.

    persistent buffer warned

    if isempty(buffer)
        buffer = strings(0, 1);
    end
    if isempty(warned)
        warned = false;
    end

    max_buffered = 5000;    % cap so a long run with a locked log cannot grow without bound
    n_attempts = 3;
    retry_pause_s = 0.05;

    try
        if nargin >= 2
            buffer(end+1, 1) = string(message); %#ok<AGROW>
        end

        if isempty(buffer)
            return
        end

        if numel(buffer) > max_buffered
            n_dropped = numel(buffer) - max_buffered;
            buffer = buffer(end-max_buffered+1:end);
            buffer(1) = sprintf("[%d earlier log lines dropped: log file stayed unavailable]", n_dropped);
        end

        % Try to append everything that is pending. Only clear the buffer once
        % the write has actually succeeded.
        for attempt = 1:n_attempts
            fid = fopen(log_path, "a");
            if fid >= 0
                try
                    fprintf(fid, "%s\n", buffer);
                    fclose(fid);
                    buffer = strings(0, 1);
                    return
                catch
                    % The handle opened but the write failed. Close it and let
                    % the next attempt or the next call try again.
                    try
                        fclose(fid);
                    catch
                    end
                end
            end
            if attempt < n_attempts
                pause(retry_pause_s);
            end
        end

        if ~warned
            warning("Could not write to %s (it may be open in another program). Log lines are buffered and will be retried.", log_path);
            warned = true;
        end
    catch
        % Never let the logger itself stop a simulation.
    end
end
