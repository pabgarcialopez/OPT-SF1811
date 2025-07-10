
% Obtain data
hw3data;

n = 8;
m = 25;
r = 3:0.25:9;

% Matrices in which results will be stored.
results_x = zeros(n, m);
results_y = zeros(n, m);
results_u = zeros(m, 1);
results_v = zeros(m, 1);
results_sigma = zeros(m, 1);
results_mu = zeros(m, 1);

Aeq = [-mu'; ones(8,1)'];

for i = 1:25
    beq = [-r(i); 1];
    [x, fval, exitflag, output, lambda] = quadprog(C, [], [], [], Aeq, beq, zeros(n, 1), Inf*ones(n, 1));

    % Write results down
    results_x(:, i) = x;
    results_y(:, i) = lambda.lower;
    results_u(i) = lambda.eqlin(1);
    results_v(i) = lambda.eqlin(2);
    results_sigma(i) = sqrt(x'*C*x);
    results_mu(i) = mu'*x;
end

% Creating a line plot
plot(results_sigma, results_mu, '-o'); 
title('Line Plot of σ vs. µ');
xlabel('σ');
ylabel('µ');






