function A_hat = optimizeSingleMatrix(constraints,objective,A,gamma,Gamma)
% Solve a single matrix optimization problems with constraints on
% singular values -- straightforward implementation of Algorithm 1.
%
% Input:
% constraints - constraint on the matrix A (yalmip)
% objective - objective in the matrix A (yalmip)
% A - the sdpvar corresponding to the matrix A (assigned to its initialization)
% gamma, Gamma - sdpvar corresponding to the bounds on singular values
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% parameters
tol = 1e-5; % stopping criteria
max_iter = 100; % maximal number of iterations
opts = sdpsettings('solver','mosek','verbose',0); % yalmip+solver options
verbose = true;

% init
iter = 0;
cont = true;
N = size(A,1);

% optimize
fprintf('----------------------------------------------------------------------------------------------------------------------\n');
while cont
    % iterate
    iter = iter+1;
    A_prev = double(A); % result of previous iteration
    
    % bound the maximal singular value (eq. (5) in the paper)
    u_bound = ([Gamma*eye(N), A; A', Gamma*eye(N)]>=0); % LMI
    
    % choose the rotation R of C_gamma (Lemma 1 in the paper)
    [U,~,V] = svd(A_prev);
    R = U*V'; % rotation of the polar decomposition of A
    RA = R'*A; % rotate A instead of rotating C_gamma (A \in R*C_gamma <==> R'*A \in C_gamma)
    
    % bound the minimal singular value (eq. (6) and (11.d) in the paper)
    l_bound = ((RA+RA')/2 >= gamma*eye(N)); % LMI
    
    % optimize
    res = solvesdp(constraints+u_bound+l_bound, objective, opts);
    A_hat = double(A);
    
    % display progress
    if verbose
        [~,S,~] = svd(A_hat);
        err = norm(A_hat-A_prev,'fro');
        fprintf('iter: %d \tobj: %7.3g \t||A-A_prev||: %7.3g \ts_min: %7.3g \ts_max: %7.3g \t\t(s_max/s_min: %7.3g) \n', iter, double(objective), err, S(end), S(1), S(1)/S(end));
    end
    
    % stopping criteria
    if (err<=tol)
        fprintf('stopping (err<=tol)\n')
        cont = false;
    elseif  (iter>=max_iter)
        fprintf('stopping: (iter>max_iter)\n')
        cont = false;
    end
end
fprintf('---------------------------------------------------------------------------------------------\n');
