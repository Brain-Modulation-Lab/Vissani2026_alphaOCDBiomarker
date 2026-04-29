function DB_Avg = consolidate_within_repetition(DB, group_keys)
% General function to average within repetitions for any set of group keys.
% Inputs:
%   - DB: a table with data
%   - group_keys: a string array of column names to group by
% Output:
%   - DB_Avg: grouped and averaged table with column `n_repetitions`

%     % Convert all keys to categorical if needed
%     for k = 1:numel(group_keys)
%         col = DB.(group_keys(k));
%         if iscellstr(col) || isstring(col)
%             DB.(group_keys(k)) = categorical(cellstr(col));
%         end
%     end

    [G, idx] = findgroups(DB(:, group_keys));
    n_groups = height(idx);
    DB_Avg = table();
    varnames = DB.Properties.VariableNames;

    for g = 1:n_groups
        group_data = DB(G == g, :);
        n_reps = height(group_data);

        % Assign group identifiers
        for k = 1:numel(group_keys)
            DB_Avg.(group_keys(k))(g,1) = idx.(group_keys(k))(g);
        end

        % Store number of repetitions
        DB_Avg.n_repetitions(g,1) = n_reps;

        for v = 1:numel(varnames)
            varname = varnames{v};
            if any(strcmp(varname, group_keys))
                continue;
            end

            col = group_data.(varname);

            try
                if isnumeric(col) || islogical(col)
                    if iscell(col)
                        col_mat = cat(1, col{:});
                    elseif ismatrix(col) && size(col, 1) == height(group_data)
                        col_mat = col;
                    elseif isvector(col)
                        DB_Avg.(varname)(g,1) = mean(col, 'omitnan');
                        continue;
                    else
                        continue;
                    end

                    avg_val = mean(col_mat, 1, 'omitnan');

                    if numel(avg_val) > 1
                        DB_Avg.(varname){g,1} = avg_val;
                    else
                        DB_Avg.(varname)(g,1) = avg_val;
                    end

                elseif iscell(col) && isnumeric(col{1})
                    col_mat = cat(1, col{:});
                    avg_val = mean(col_mat, 1, 'omitnan');
                    DB_Avg.(varname){g,1} = avg_val;

                elseif iscellstr(col) || isstring(col) || iscategorical(col) || ischar(col)
                    % Convert char to string if needed
                    if ischar(col)
                        col = string({col});
                    elseif iscellstr(col)
                        col = string(col);
                    elseif iscategorical(col)
                        col = string(cellstr(col));
                    end

                    % Check uniformity
                    if all(col == col(1))
                        DB_Avg.(varname)(g,1) = col(1);
                    else
                        DB_Avg.(varname)(g,1) = missing;
                    end
                end
            catch
                continue;
            end
        end
    end



% Convert cell columns with uniform row length into matrices
cell_vars = varfun(@iscell, DB_Avg, 'OutputFormat', 'uniform');
for v = find(cell_vars)
    varname = DB_Avg.Properties.VariableNames{v};
    try
        DB_Avg.(varname) = cellcol2mat(DB_Avg.(varname));
    catch
        % If conversion fails (e.g. unequal lengths), leave as cell
    end
end

end



function mat = cellcol2mat(cellcol)
% Convert a column of cells (each containing a row vector) into a matrix

    if ~iscell(cellcol)
        error('Input must be a cell array.');
    end
    if isempty(cellcol)
        mat = [];
        return;
    end
    % Check that all cells have the same length
    lengths = cellfun(@length, cellcol);
    if any(lengths ~= lengths(1))
        error('Not all rows in the cell array are the same length.');
    end

    mat = vertcat(cellcol{:});
end
