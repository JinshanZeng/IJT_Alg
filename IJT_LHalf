%@ Edited by Jinshan Zeng, School of Mathematics&Stastics, Xi'an Jiaotong University, 
% email: jsh.zeng@gmail.com, date: Sep. 27, 2014
% Iterative Jumping Thresholding algorithm for L_{1/2} regularization problem:
% min ||A*x - y||^2 + lambda ||x||_{0.5}^{0.5}
function [x, cpu_time] = IJT_LHalf(A, y, lambda, x0, MaxIterNum)
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
    id = find(abs(b)>(54)^(1/3)/4*lambda_mu^(2/3));
    x = zeros(n,1);
    x(id) = real(2/3*b(id).*(1 + cos(2/3*(pi - acos(lambda_mu/8*((abs(b(id))/3).^(-1.5)))))));
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
