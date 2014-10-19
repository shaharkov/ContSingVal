%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script:
% Simple volumetric mesh optimization problems with constraints on singular values.
% An application of Algorithm 2 to a few example optimization problems,
% generating some of the example deformation of Figure 2 in the paper.
% See section 6.1 for additional details.
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
[X, tri, anchors, anchor_coords] = generateBarExample;

%% parameters
C = 2; % bound on conformal distortion

%% AAAP unconstrained
s{1} = Solver(Problem(X,tri,0,ObjectiveEnum.AAAP,SpaceEnum.NONE)); % setup problem
s{1}.problem.auxConstraints = (s{1}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{1}.problem.setFrames([]); % initialize rotations
s{1}.solve; % solve
res.aaap_unconstrained.Y = double(s{1}.problem.Y); % store result

%% AAAP BD
s{2} = Solver(Problem(X,tri,C,ObjectiveEnum.AAAP,SpaceEnum.BD)); % setup problem
s{2}.problem.auxConstraints = (s{2}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{2}.problem.setFrames(s{1}.problem); % initialize rotations (with the result of unconstrained AAAP)
s{2}.solve;
res.aaap_bd.Y = double(s{2}.problem.Y);

%% AAAP BSI
s{3} = Solver(Problem(X,tri,sqrt(C),ObjectiveEnum.AAAP,SpaceEnum.BSI)); % setup problem
s{3}.problem.auxConstraints = (s{3}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{3}.problem.setFrames(s{1}.problem); % initialize rotations (with the result of unconstrained AAAP)
s{3}.solve;
res.aaap_bsi.Y = double(s{3}.problem.Y);

%% ARAP unconstrained
s{4} = Solver(Problem(X,tri,0,ObjectiveEnum.ARAP,SpaceEnum.NONE)); % setup problem
s{4}.problem.auxConstraints = (s{4}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{4}.problem.setFrames(s{1}.problem); % initialize rotations (with the result of AAAP)
s{4}.solve;
res.arap_unconstrained.Y = double(s{4}.problem.Y);

%% ARAP BD
s{5} = Solver(Problem(X,tri,C,ObjectiveEnum.ARAP,SpaceEnum.BD)); % setup problem
s{5}.problem.auxConstraints = (s{5}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{5}.problem.setFrames(s{4}.problem); % initialize rotations (with the result of unconstrained ARAP)
s{5}.solve;
res.arap_bd.Y = double(s{5}.problem.Y);

%% ARAP BSI
s{6} = Solver(Problem(X,tri,sqrt(C),ObjectiveEnum.ARAP,SpaceEnum.BSI)); % setup problem
s{6}.problem.auxConstraints = (s{6}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{6}.problem.setFrames(s{4}.problem); % initialize rotations (with the result of unconstrained ARAP)
s{6}.solve;
res.arap_bsi.Y = double(s{6}.problem.Y);

%% EQC (for additional details see example_ExtremalQuasiConformal.m)
initC=50;
tolC=1e-2;
s{7} = Solver(Problem(X,tri,0,ObjectiveEnum.FEASIBILITY,SpaceEnum.BD)); % setup problem
s{7}.problem.auxConstraints = (s{7}.problem.Y(anchors,:)==anchor_coords); % set positional constraints
s{7}.problem.setFrames([]);
curC=initC;
lastC=inf;
while (lastC-curC)>tolC,
    s{7}.problem.C=curC;
    s{7}.solve();
    % update C
    lastC=curC;   
    curC=max(s{7}.problem.distortion);
    assert(curC<=lastC,'C is not monotone');
    fprintf(2,'Updating C=%g\n',curC);
end
fprintf(2,'Done. C=%g\n',curC);
res.eqc.Y = double(s{7}.problem.Y);

%% plot everything
boundary_tri = getBoundaryFaces(tri); % triangulation of the volume's boundary

figure;
fields = fieldnames(res);
for i = 1:length(fields)
    exp = fields{i};
    subplot(2,4,i);
    patch('faces',boundary_tri,'vertices',res.(exp).Y,'facecolor','c')
    cameratoolbar;
    cameratoolbar('SetCoordSys','y');
    view(-10,20)
    axis equal;
    axis off;
    title(exp);
    set(get(gca,'title'),'Interpreter','none')
end