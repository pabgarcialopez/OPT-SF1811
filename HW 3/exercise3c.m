

% Obtain data
hw3data;

n = 8;
m = 25;
r = 3:0.25:9;

% Matrices in which results will be stored.
results_x_1 = zeros(n, m);
results_y_1 = zeros(n, m);
results_u_1 = zeros(m, 1);
results_v_1 = zeros(m, 1);
results_sigma_1 = zeros(m, 1);
results_mu_1 = zeros(m, 1);

Aeq = [-mu'; ones(8,1)'];

for i = 1:25
    beq = [-r(i); 1];
    [x, fval, exitflag, output, lambda] = quadprog(C, [], [], [], Aeq, beq, zeros(n, 1), Inf*ones(n, 1));

    % Write results down
    results_x_1(:, i) = x;
    results_y_1(:, i) = lambda.lower;
    results_u_1(i) = lambda.eqlin(1);
    results_v_1(i) = lambda.eqlin(2);
    results_sigma_1(i) = sqrt(x'*C*x);
    results_mu_1(i) = mu'*x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve modified problem

n = 8;
m = 25;
r = 3:0.25:9;
results_x_2 = zeros(n, m);
results_y_2 = zeros(n, m);
results_u_2 = zeros(m, 1);
results_v_2 = zeros(m, 1);
results_sigma_2 = zeros(m, 1);
results_mu_2 = zeros(m, 1);

A = -mu';
Aeq = ones(1, n);

for i = 1:25
    [x, fval, exitflag, output, lambda] = quadprog(C, [], A, -r(i), Aeq, 1, zeros(n, 1), Inf*ones(n, 1));

    % Write results down
    results_x_2(:, i) = x;
    results_y_2(:, i) = lambda.lower;
    results_u_2(i) = lambda.ineqlin;
    results_v_2(i) = lambda.eqlin(1);
    results_sigma_2(i) = sqrt(x'*C*x);
    results_mu_2(i) = mu'*x;
end


% Creating a line plot
plot(results_sigma_1, results_mu_1, '-o'); 
hold on;

% Creating the second line plot
plot(results_sigma_2, results_mu_2, '-s');

xlabel('σ');
ylabel('µ');
title('Line Plot of σ vs. µ');

% Adding a legend
legend('Original problem', 'Modified problem');

hold off;





