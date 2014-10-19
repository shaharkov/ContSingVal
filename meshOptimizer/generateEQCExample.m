function [X, tri, anchors, anchor_coords] = generateEQCExample(N)
% generates the volumetric mesh and constraints of Figure 1 in the paper
% (extremal quasiconformal mapping).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate vertices
[x,y,z]=meshgrid(1:N,1:N,1:N);
x=x(:);y=y(:);z=z(:);
X=[x y z];

% important coordinates
ubound = [N N N];
ctr = (ubound+1)/2;
ctr_round = floor(ctr);

% generate mesh
tri = delaunay(X);

% constraints
anchors = find((x==ubound(1) | x==1) & (y==ubound(2) | y==1) & (z==ubound(3) | z==1)); % corners
anchor_coords = X(anchors,:);
anchor_coords = anchor_coords + 0.2*N*randn(size(anchor_coords));
