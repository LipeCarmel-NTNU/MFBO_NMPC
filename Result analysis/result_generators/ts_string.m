function ts = ts_string(dt)
%TS_STRING Filename timestamp string yyyyMMdd_HHmmss from a datetime column.
ts = string(datetime(dt, 'Format', 'yyyyMMdd_HHmmss'));
end
