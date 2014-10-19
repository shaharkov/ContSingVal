function [X, tri, anchors, anchor_coords] = generateBarExample
% generates the volumetric mesh and constraints of Figure 2 in the paper.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate vertices
Nb = 5;Nh = 30;
[x,y,z]=meshgrid(1:Nh,1:Nb,1:Nb);
x=x(:);y=y(:);z=z(:);
X=[x y z];

% important coordinates
ubound = [Nh Nb Nb];
ctr = (ubound+1)/2;
ctr_round = floor(ctr);

% generate mesh
tri = delaunay(X);

% constraints
anchors_fix = find(x<=2);
anchors_fix_coords = X(anchors_fix,:);
anchors_move = find(x>=ubound(1));
R = getRotation3D(70*pi/180,45*pi/180,0);
anchors_move_coords = bsxfun(@plus,bsxfun(@plus,bsxfun(@minus,X(anchors_move,:),ctr_round)*R*1,ctr_round),ubound.*[-1/2 0 -1/2]);
anchors = [anchors_fix;anchors_move];
anchor_coords = [anchors_fix_coords;anchors_move_coords];