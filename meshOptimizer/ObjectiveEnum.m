classdef ObjectiveEnum<handle
% Objective type
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

    enumeration
        ARAP; % As-Rigid-As-Possible
        AAAP; % As-Affine-As-Possible
        DIRICHLET; % Dirichlet energy
        FEASIBILITY; % feasibility related functional (used for computing extremal quasiconformal mappings)
    end   
end

