classdef Stopper
% Helper class for measuring times.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        t_start;
        t_end;
    end   
    methods
        function obj=Stopper(varargin)
            obj.t_start = tic;
            if length(varargin)>=1,
                fprintf(varargin{:});
            end;
        end
        function t_end=stop(obj,varargin)
            t_end = toc(obj.t_start);
            obj.t_end = t_end;
            if length(varargin)==0,
                fprintf('(%.3g secs)\n',t_end);
            end;
        end
    end
end

