% Define inversions and meshing options
options         = optimset('Display', 'off');
options2        = optimset('MaxIter', 1000);
options3.dhmax = 10;
% Define other variable for jRi, inversion, and terminations
covd2       = repmat(diag(ones(1,np)),2,2);
lambdas      = 10.^([linspace(-2,2,15)]);
tol         = 0.1;                  % Termination tolerance