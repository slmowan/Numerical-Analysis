% create original sys
n = 6;
v = (1:n) .^ 2;
A = vander(v);

% create E (perturb A by abs(A))
E = abs(A);
% create b
b = A * ones(n, 1);
% create e (perturb b by abs(b))
e = abs(b);

% create the perturbed system
tau = 10^(-10);
A_perturbed = A + tau*randn(n,n) .* E;
b_perturbed = b + tau*randn(n,1) .* e;


% solve the system
% forward err first
y = A_perturbed \ b_perturbed; % \ is inverse, and this is forward error
% calculate forward errro
x = ones(n,1); % exact solution
forward_err = norm(y-x)/norm(x);
disp("forward error is :")
disp(forward_err)

% cauculate the bound, RHS of theorem 2, because we want to compare the cound with forward err to analyze
% calculate numerator first
num = tau*cond(A) * (norm(E)/norm(A) + norm(e)/norm(b));
disp('numerator is :')
disp(num);
% then calculate denomenator
den = 1 - tau*norm(inv(A)) * norm(E);
disp('denomenator is :')
disp(den);
% finally, calculate bound
bound = num/den;
disp('bound is :')
disp(bound);

% now, backward err, r = b - A*y
r = b - A*y;
backward_err = norm(r)/(norm(E)*norm(y) + norm(e));
disp('backward error is :')
disp(backward_err);

% forward bound (provided by theorem 2)