function ts = timestamp_compact()
%TIMESTAMP_COMPACT Compact yyyymmdd_HHMMSS timestamp in Europe/Oslo time.
    t = datetime("now", "TimeZone", 'Europe/Oslo');
    ts = char(datestr(t, "yyyymmdd_HHMMSS"));
end
