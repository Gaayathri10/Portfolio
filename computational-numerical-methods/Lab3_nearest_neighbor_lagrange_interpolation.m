clc; clear; close all;

% (a) Find f(x) values for x = 5.6, 9.2, and 13.05 using nearest neighbor interpolation

% the data set
x = [2, 4, 7, 10, 12, 15];
f = [5, 18, 200, 310, 320, 450];

x_query = [5.6, 9.2, 13.05];

% Using knnsearch to find the nearest neighbor indices
Idx = knnsearch(x', x_query');

% Fetch corresponding f(x) values
f_nn = f(Idx);

disp('Nearest neighbor interpolation results:');
disp(table(x_query', f_nn', 'VariableNames', {'x_query', 'f_nn'}));



% (b) Implement a MATLAB function for Lagrange interpolation with k nearest points

function result = lagrange_interpolation(x, f, x_query, k)
    % Finding the k closest points
    [~, sorted_idx] = sort(abs(x - x_query));  % Sort by distance
    nearest_x = x(sorted_idx(1:k));  % Select k nearest x values
    nearest_f = f(sorted_idx(1:k));  % Select corresponding f(x) values

    result = 0;

    % Computing the Lagrange interpolation
    for i = 1:k
        term = nearest_f(i);
        for j = 1:k
            if j ~= i
                term = term * (x_query - nearest_x(j)) / (nearest_x(i) - nearest_x(j));
            end
        end
        result = result + term;
    end
end

% (c) Running the Lagrange interpolation for x = 5.6, 9.2, and 13.05 using k = 2, 3, and 4

k_values = [2, 3, 4];

disp('Lagrange Interpolation Results:');
for i = 1:length(x_query)
    for j = 1:length(k_values)
        result = lagrange_interpolation(x, f, x_query(i), k_values(j));
        fprintf('Interpolated f(%.2f) using %d nearest points: %.4f\n', x_query(i), k_values(j), result);
    end
end
