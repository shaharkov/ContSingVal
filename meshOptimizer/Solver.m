classdef Solver<handle
% Wrapper class for mesh optimization problem with constraints
% on singular values. Implements Algorithm 2 in the paper.
% Can be used for various mesh optimization problems
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
        problem; % problem instance       
        maxIter = 50; % maximal number of optimization iterations
        tolY = 1e-2; % stopping criteria w.r.t. vertices
        tolObj = 1e-2; % stopping criteria w.r.t. objective
		feasibilityStopOnFeasibility=true; % whether to stop whenever feasible (just with feasibility functional) 
        plotLogPerIter=true; % whether to plot graphs every iteration
		lastY; % last configuration of vertices
        lastObj; % last objective value
        log=[]; % optimization log
    end
    
    methods
        function obj=Solver(problem)
			% setup solver
            obj.problem=problem;
        end
        
        function solve(obj)
			% main function -- implementing Algorithm 2.
            fprintf('\n====================================================\n');
            % visualize initial state
            figure(200);
            obj.problem.visualize();
            drawnow; pause(0.01);
            
            % log initial state
            obj.problem.dispReport();
            obj.updateLog(0);
            
            % initialize solver loop
            obj.lastY=obj.problem.X;
            lastObj=inf;
            
            % iterate
            iter=1;
            while(iter<=obj.maxIter)
                t_start = tic;
                fprintf('\n==========================\n');
                fprintf('iter %d\n',iter);
                
                % main optimization block
                obj.problem.optimize();
                obj.problem.updateDistortion();
                obj.problem.dispReport();
                obj.updateLog(iter);
                
                % visualize progress
                figure(200); cla;
                obj.problem.visualize();
                if obj.plotLogPerIter,
                    figure(201);
                    obj.plotLog();
                end;
                drawnow; pause(0.01);
                
                % update rotations
                Y=double(obj.problem.Y);
                obj.problem.updateFrames();
                
                % conclude iteration
                t_stop = toc(t_start);
                fprintf('iteration time: %.2f secs\n', t_stop);
                fprintf('==========================\n');
                iter=iter+1;
                
                % stopping criteria
                if norm(Y(:)-obj.lastY(:),'inf')<obj.tolY
                    fprintf('Stopping (tolY)\n');
                    break;
                end
                if abs((obj.lastObj-obj.problem.objVal)/obj.problem.objVal)<obj.tolObj
                    fprintf('Stopping (tolObj)\n');
                    break;
                end
                if (obj.problem.objType==ObjectiveEnum.FEASIBILITY) && obj.feasibilityStopOnFeasibility
                    fprintf('Stopping (Feasible)\n');
                    break;
                end
                
                % save previous iteration
                obj.lastY=Y;
                obj.lastObj=obj.problem.objVal;
            end
            fprintf('Done.\n');
        end
        
		%% functions for logging and plotting
        function updateLog(obj,iter)
            indLast=length(obj.log);
            ind=indLast+1;
            obj.log(ind).iter=iter;
            obj.log(ind).time=now;
            obj.log(ind).Y=double(obj.problem.Y);
            obj.log(ind).F=cat(3,obj.problem.faces.F);
            obj.log(ind).C=obj.problem.C;
            obj.log(ind).objVal=obj.problem.objVal;
            obj.log(ind).distortion=obj.problem.distortion;
            obj.log(ind).flipped=obj.problem.flipped;
            obj.log(ind).problem_log=obj.problem.log;
            fprintf('Updating log item #%d (%s)\n',ind,datestr(now));
        end
        function restoreFromLog(obj,ind)
            assign(obj.problem.Y,obj.log(ind).Y);
            for i=1:size(obj.log(ind).F,3),
                obj.problem.faces(i).F=obj.log(ind).F(:,:,i);
            end
            obj.problem.C=obj.log(ind).C;
            obj.problem.objVal=obj.log(ind).objVal;
            obj.problem.distortion=obj.log(ind).distortion;
            obj.problem.flipped=obj.log(ind).flipped;
            obj.problem.active=obj.log(ind).active;
            obj.problem.log=obj.log(ind).problem_log;
        end
        function plotLog(obj,range)
            clf;
            % prepare data
            iter=cat(1,obj.log.iter);
            objVal=cat(1,obj.log.objVal);
            for i=1:length(obj.log),
                maxDist(i)=max(obj.log(i).distortion);
                numFlips(i)=nnz(obj.log(i).flipped);
            end
            counter=1:length(iter);
            if ~exist('range','var')
                range=counter(1:end);
            end
            
            % plot
            AX = gca;
            plot(counter(range),objVal(range));
            addaxis(counter(range),maxDist(range));
            addaxislabel(1,'Objective');
            addaxislabel(2,'Maximal distortion');
            if any(numFlips)
                addaxis(counter(range),numFlips(range),[min(numFlips) max(numFlips)+1]);
                addaxislabel(3,'Number of flips');
            end
            
            % add vertical lines for resets
            ind=range(iter(range)==0);
            for i=1:length(ind),
                line([ind(i) ind(i)],get(AX(1),'YLim'),'Color',[.5 .5 .5],'LineStyle',':')
            end
        end  
    end
end

