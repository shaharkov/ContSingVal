%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script:
% Implementing an algorithm for computing extremal quasiconformal mappings 
% of volumetric meshes (i.e., minimizing maximal conformal distortion).
% The code reproduces the example presented in Figure 1. 
% See section 6.1 for additional details.
% Modify N to change the cube size.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
rng(1)
close all
clear
initialize;
yalmip('clear')

%% prepare example
N = 15; % set cube size (NxNxN)
[X, tri, anchors, anchor_coords] = generateEQCExample(N);

%% setup problem
initC=50; % initial bound on conformal distortion
tolC=1e-2; % stopping criteria for conformal distortion convergence
s = Solver(Problem(X,tri,0,ObjectiveEnum.FEASIBILITY,SpaceEnum.BD)); % setup problem
s.problem.auxConstraints = (s.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s.problem.setFrames([]); % initialize rotations (identity)

%% optimize
curC=initC;
lastC=inf;
while (lastC-curC)>tolC,
    % optimize with fixed C
    s.problem.C=curC;
    s.solve();
    % update C
    lastC=curC;
    curC=max(s.problem.distortion);
    assert(curC<=lastC,'C is not monotone');
    fprintf(2,'Updating C=%g\n',curC);
end
fprintf(2,'Done. C=%g\n',curC);
res.eqc.Y = double(s.problem.Y);

%% plot everything
boundary_tri = getBoundaryFaces(tri); % triangulation of the volume's boundary

figure;
patch('faces',boundary_tri,'vertices',res.eqc.Y,'facecolor','c')
cameratoolbar;
cameratoolbar('SetCoordSys','y');
view(-10,20)
axis equal;
axis off;
title('extremal quasiconformal');
set(get(gca,'title'),'Interpreter','none')
