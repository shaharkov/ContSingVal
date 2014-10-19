function initialize
% Initialization -- add required toolboxes to the path.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YalmipPath = '../yalmip/'; % <-------- update with your YALMIP path
MosekPapth = 'C:/Program Files/Mosek/7/toolbox/r2013a'; % <-------- update with your MOSEK path
% set path
if isempty(whos('global','path_def'))
    fprintf('- Adding toolbox paths\n');
    addpath('meshOptimizer');
    addpath('addaxis'); % addaxis toolbox is used for plotting progress
    addpath(genpath(YalmipPath));
    addpath(genpath(MosekPapth));
    global path_def
end
