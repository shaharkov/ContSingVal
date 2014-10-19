%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script:
% Simple single matrix optimization problems with constraints on singular values.
% An application of Algorithm 1 to a few example optimization problems,
% aiming to illustrate the implementation and usage of the theory presented
% in the paper (section 4 and 5).
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% init
rng(1)
close all
clear
initialize
yalmip('clear')


%% setup variables
N = 5; % matrix dimension
A = sdpvar(N,N,'full');
gamma = sdpvar;
Gamma = sdpvar;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find closest matrix to a random matrix with SVs bounded in [0.9,1.1]
% (see problem presented in eq. (1) in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
objective = norm(A-randn(N), 'fro'); % ||A-rand||_fro
constraints = (0.9<=gamma)+(Gamma<=1.1);
assign(A,eye(N)); % init A to be identity
A_hat = optimizeSingleMatrix(constraints,objective,A,gamma,Gamma)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find closest non-negative matrix to a random matrix with condition number
% (sigma_max/sigma_min) bounded by 1.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
objective = norm(A-randn(N), 'fro'); % ||A-rand||_fro
constraints = (Gamma<=1.5*gamma);
assign(A,eye(N)); % init A to be identity
A_hat = optimizeSingleMatrix(constraints,objective,A,gamma,Gamma)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a matrix with a given random first row that minimizes that spread
% of singular values (sigma_max-sigma_min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
objective = (Gamma-gamma);
constraints = (A(1,:)==randn(1,N));
assign(A,eye(N)); % init A to be identity
A_hat = optimizeSingleMatrix(constraints,objective,A,gamma,Gamma)
