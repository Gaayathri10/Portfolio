
function romberg_integration()
    % Defining the function f(t) = Q(t) * c(t)
    f = @(t) (2 + 4*(cos(0.2*t)).^2) .* (5*exp(-0.4*t) + 2*exp(0.13*t));
    
    % Defining the time intervals
    intervals = [1 7; 3 9; 2 6];
    
    % Tolerance of 0.1%
    tol = 0.001;
    
    % Solve for each interval
    for i = 1:size(intervals, 1)
        a = intervals(i, 1);
        b = intervals(i, 2);
        
        fprintf('\nCalculating for interval t1 = %g to t2 = %g minutes\n:', a, b);
        
        % Calling romberg function with adaptive stopping
        [R, n] = romberg_adaptive(f, a, b, tol);
        fprintf('Final result: M = %.6f mg\n', R(n, n));
    end
end

function [R, n] = romberg_adaptive(f, a, b, tol)
    % Initialize variables
    max_iter = 20; % Maximum no of iterations to prevent infinite loops
    R = zeros(max_iter, max_iter);
    n = 1;
    h = b - a; % step size
    
    % First trapezoidal rule
    R(1,1) = h * (f(a) + f(b)) / 2;
    
    % Romberg iteration
    for n = 2:max_iter
        % Composite trapezoidal rule 
        h = h / 2;
        subtotal = 0;
        for i = 1:2^(n-2)
            subtotal = subtotal + f(a + (2*i-1)*h);
        end
        R(n,1) = R(n-1,1)/2 + h*subtotal;
        
        % Richardson extrapolation
        for m = 2:n
            R(n,m) = (4^(m-1)*R(n,m-1) - R(n-1,m-1)) / (4^(m-1) - 1);
        end
        
        % Check for convergence
        if n > 1 && abs(R(n,n) - R(n-1,n-1)) < tol * abs(R(n,n))
            break;
        end
    end
    
    % Display the Romberg table with 6 decimal places
    fprintf('\nRomberg Table:\n');
    h_values = (b-a) ./ (2.^(0:n-1))';  % Step sizes
    
    % Creating a cell array 
    col_names = ["n", "h"];
    for j = 1:n
        col_names = [col_names, "O(h^" + num2str(2*j) + ")"];
    end
    
    % Creating a matrix of the romberg table 
    R_table = [ (0:n-1)', h_values, R(1:n, 1:n) ];
    
    % Formatting the values to be 6 decimal values 
    R_table = arrayfun(@(x) sprintf('%.6f', x), R_table, 'UniformOutput', false);
    
    % Converting the table and display
    R_table = cell2table(R_table, 'VariableNames', col_names);
    disp(R_table);
end
