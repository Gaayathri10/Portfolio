A = [3, -1, -1; 
     6, 7, 2; 
    -5, 1, 8];

b = [1; 25; 35];

rank_A = rank(A); % Rank number
cond_A = cond(A); % Condition number

% Displaying results for Q1 w/15 decimal places
disp(['Rank of A: ', num2str(rank_A, 15)]);
disp(['Condition number of A: ', num2str(cond_A, 15)]);

% Q2: Solving using matrix inversion (if justified)
if cond_A < 1e10  % Checking if inversion is numerically stable
    x_inverse = inv(A) * b;  % Using matrix inversion
else
    disp('Matrix is ill-conditioned. Inversion may not be reliable.');
    x_inverse = NaN(3,1); % Avoid undefined variable error
end

% Solving using MATLAB backslash operator
x_backslash = A \ b;

% Display results for Q2 w/15 decimals
disp('Solution using matrix inversion:');
disp(['x = ', num2str(x_inverse(1), 15)]);
disp(['y = ', num2str(x_inverse(2), 15)]);
disp(['z = ', num2str(x_inverse(3), 15)]);

disp('Solution using backslash operator:');
disp(['x = ', num2str(x_backslash(1), 15)]);
disp(['y = ', num2str(x_backslash(2), 15)]);
disp(['z = ', num2str(x_backslash(3), 15)]);

% Q3: LU decomposition (Doolittle)
n = size(A,1);
L = eye(n); % Lower triangular matrix initialized as identity
U = zeros(n); % Upper triangular matrix initialized as zeros

% LU decomposition
for i = 1:n
    for k = i:n
        U(i,k) = A(i,k) - L(i,1:i-1) * U(1:i-1,k);
    end
    for k = i+1:n
        L(k,i) = (A(k,i) - L(k,1:i-1) * U(1:i-1,i)) / U(i,i);
    end
end

% Forward substitution to solve L*y = b
y = zeros(n,1);
for i = 1:n
    y(i) = (b(i) - L(i,1:i-1) * y(1:i-1)) / L(i,i);
end

% Backward substitution to solve U*x = y
x_LU = zeros(n,1);
for i = n:-1:1
    x_LU(i) = (y(i) - U(i,i+1:n) * x_LU(i+1:n)) / U(i,i);
end

% Displaying results
disp('Solution using LU decomposition (Doolittle):');
disp(['x = ', num2str(x_LU(1), 15)]);
disp(['y = ', num2str(x_LU(2), 15)]);
disp(['z = ', num2str(x_LU(3), 15)]);

%Q4 Solving using Jacobi method 
x_jacobi = zeros(3,1); % Initial guess
tol = 0.0001; % Convergence tolerance
max_iter = 100; % Maximum iterations

for iter = 1:max_iter
    x_new = zeros(3,1);
    x_new(1) = (1 + x_jacobi(2) + x_jacobi(3)) / 3;
    x_new(2) = (25 - 6*x_jacobi(1) - 2*x_jacobi(3)) / 7;
    x_new(3) = (35 + 5*x_jacobi(1) - x_jacobi(2)) / 8;
    
    % Computing relative error
    error = max(abs((x_new - x_jacobi) ./ x_new));
    x_jacobi = x_new;
    
    if error < tol
        break;
    end
end

% Displaying results for Q4
disp('Solution using Jacobi method:');
disp(['x = ', num2str(x_jacobi(1), 15)]);
disp(['y = ', num2str(x_jacobi(2), 15)]);
disp(['z = ', num2str(x_jacobi(3), 15)]);
disp(['Iterations: ', num2str(iter)]);


%Q5 Gauss-sidel Method 
x_gs = zeros(3,1); % Initial guess 
tol = 0.0001; % Convergence tolerance
max_iter = 100; % Maximum iterations

for iter = 1:max_iter
    x_old = x_gs; % Storing previous values for error calc
    
    % Updating each variable using Gauss-Seidel method
    x_gs(1) = (1 + x_gs(2) + x_gs(3)) / 3;
    x_gs(2) = (25 - 6*x_gs(1) - 2*x_gs(3)) / 7;
    x_gs(3) = (35 + 5*x_gs(1) - x_gs(2)) / 8;
    
    % Computing relative error
    error = max(abs((x_gs - x_old) ./ x_gs));
    
    if error < tol
        break;
    end
end

% Displaying results for Gauss-Seidel method
disp('Solution using Gauss-Seidel method:');
disp(['x = ', num2str(x_gs(1), 15)]);
disp(['y = ', num2str(x_gs(2), 15)]);
disp(['z = ', num2str(x_gs(3), 15)]);
disp(['Iterations: ', num2str(iter)]);
