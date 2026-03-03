function bytes = get_precision_bytes(precision)
%% Return the number of bytes per element for a given precision
switch lower(precision)
    case {'float32', 'single'}, bytes = 4;
    case {'float64', 'double'}, bytes = 8;
    case {'int32'}, bytes = 4;
    case {'int16'}, bytes = 2;
    case {'uint8'}, bytes = 1;
    otherwise
        error('Unsupported precision: %s', precision);
end
end