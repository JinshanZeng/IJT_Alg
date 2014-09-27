%@ Edited by Jinshan Zeng, School of Mathematics&Stastics, Xi'an Jiaotong University, 
% email: jsh.zeng@gmail.com, date: Sep. 27, 2014
% Iterative Jumping Thresholding algorithm for L_{2/3} regularization problem:
% min ||A*x - y||^2 + lambda ||x||_{2/3}^{2/3}
function [x, cpu_time] = IJT_L2rds(A, y, lambda, x0, MaxIterNum)
%% Input %%%%%%%%%%
% A             -the measurement matrix
% y             -the measurement
% x0            -Initial Value
% lambda        - Regularization parameter
% MaxIterNum    - Maximal iteration number
%% Output %%%%%%%%%%
% x            - Recovery
% cpu_time     - running time of the algorithm
%
warning off all;
%% Start running
start_time = cputime;
n = size(A,2);
iter = 1;
r = y;
mu = 0.99*norm(A)^(-2);
% mu = 0.9*normest(A)^(-2); % Parameter algorithm 'mu'
stop = 1;
while stop==1
    u = A'*r;
    b = x0 + mu*u;
    lambda_mu = lambda*mu;  
    %% Update the solution
    id =find(abs(b)>(2/3)*(3*lambda_mu^3)^(1/4));
    x = zeros(n,1);
    a = zeros(n,1);
    a(id) = 2/sqrt(3)*(lambda_mu).^(0.25).*(cosh(acosh(27/16*b(id).^2*(lambda_mu).^(-3/2))./3)).^(0.5);
    x(id) = sign(b(id)).*real(((abs(a(id))+sqrt(2*abs(b(id))./abs(a(id))-abs(a(id)).^2))/2).^3);
    if norm(x - x0)/norm(x)>(1E-6) && iter<MaxIterNum
        stop = 1;
        %% Update the previous solution
        iter = iter + 1;
        r = y - A*x;
        x0 = x;
    else
        stop = 0;
    end
end
end_time = cputime;
cpu_time = end_time-start_time;
