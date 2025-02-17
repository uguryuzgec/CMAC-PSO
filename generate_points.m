function points = generate_points(N, dim, bounds)
    % Generate points in a specified dimensional space within given bounds
    %
    % N - Number of points to generate
    % dim - Dimensionality of the space
    % bounds - 2xD matrix with min and max bounds for each dimension

    % Calculate the number of points per dimension
    num_points_per_dim = ceil(N^(1/dim));

    % Initialize cell array to store grid vectors
    grid_vectors = cell(1, dim);
    for i = 1:dim
        linspace_vals = linspace(bounds(1, i), bounds(2, i), num_points_per_dim + 2);
        grid_vectors{i} = linspace_vals(2:end-1); % Exclude boundaries
    end

    % Preallocate array to store points
    points = zeros(N, dim);

    % Generate points using ndgrid and permute to avoid large memory consumption
    for i = 1:N
        % Convert linear index to subscript
        subs = cell(1, dim);
        [subs{:}] = ind2sub(repmat(num_points_per_dim, 1, dim), i);
        for j = 1:dim
            points(i, j) = grid_vectors{j}(subs{j});
        end
    end
end