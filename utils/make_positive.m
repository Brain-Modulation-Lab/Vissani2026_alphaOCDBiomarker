function x = make_positive(x)

    x = x(:);

    % Replace bad values
    x(~isfinite(x)) = NaN;

    % Shift to positive range
    x = x - min(x, [], 'omitnan');

    % Add small positive offset
    x = x + 0.01;

end