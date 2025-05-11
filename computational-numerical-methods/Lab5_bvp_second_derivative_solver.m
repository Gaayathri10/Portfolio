
function [x, y] = BVP2ndDriv(a, b, Ya, Db, n, pOFx, qOFx, rOFx)
    % Defineing the step size and grid
    h = (b - a)/n;
    x = linspace(a, b, n+1)'; % n+1 points

    % Initialize matrix A and RHS vector B
    A = zeros(n+1, n+1);
    B = zeros(n+1, 1);

    % Boundary condition at x = a
    A(1,1) = 1;
    B(1) = Ya;

    % Interior points in central difference
    for i = 2:n
        xi = x(i);
        pi = pOFx(xi);
        qi = qOFx(xi);
        ri = rOFx(xi);

        A(i, i-1) = 1/h^2 - pi/(2*h);
        A(i, i)   = -2/h^2 + qi;
        A(i, i+1) = 1/h^2 + pi/(2*h);
        B(i) = ri;
    end

    % Boundary condition at x = b 
    A(n+1, n) = -1/h;
    A(n+1, n+1) = 1/h;
    B(n+1) = Db;

    % Solving the system
    y = A \ B;
end

% Defining the given problem parameters
a = -1;
b = 3;
Ya = 2;
Db = -7;
n = 30;

% Defining the functions p(x), q(x), r(x)
pOFx = @(x) -2 ./ x;
qOFx = @(x) 0 .* x;
rOFx = @(x) 21 + 0 .* x;

% Solving using the defined function
[x, y] = BVP2ndDriv(a, b, Ya, Db, n, pOFx, qOFx, rOFx);

% Ploting the solution
plot(x, y, '-o', 'LineWidth', 1.5);
xlabel('x');
ylabel('y(x)');
title('Solution to the BVP');
grid on;
