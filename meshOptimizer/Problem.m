classdef Problem<handle
% A basic problem class for mesh optimization problem with constraints
% on singular values. Implements the code required for running a single
% iteration of Algorithm 2.
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
        X; % input vertices
		tris; % triangulation
		Y; % output vertices (yalmip optimization variables)
        INIT_Y; % initialization for Y
		vols; % vector of element volumes
        distortion=nan; % vector of element distortions
        flipped=nan; % vector of element flips
        faces; % class array of elements (faces)
        C; % bound on distortion (isometric, conformal, etc.)
        objType; % objective for optimization (see ObjectiveEnum)
        spaceType; % objective for optimization (see SpaceEnum)
        tau; % auxiliary variable for feasibility problems
        aux_k; % auxiliary variable for bounded scaled-isometry problems
        aux_K; % auxiliary variable for bounded scaled-isometry problems
        source_dim; % source mesh dimension
        target_dim; % source mesh dimension
        objVal=nan; % objective value
        objective; % objective yalmip optimization variable)
        objectiveConstraints; % objective related constraints (in epigraph form)
		auxConstraints; % auxiliary constraints
        recalc_objective=true; % whether objective should be recalculated each iteration (e.g., in non-linear ARAP)
        yalmip_prm; % yalmip's parameters
        log; % log
    end
    
    methods
        function obj=Problem(X,tris,C,objType,spaceType,Y)
            % initialize and setup a problem
            stopper=Stopper('Constructing problem... ');
            obj.faces=Face();
            obj.vols=computeVolumes(tris,X);
            obj.auxConstraints=[];
            obj.tau=sdpvar();
            obj.aux_k=sdpvar();
            obj.aux_K=sdpvar();
            assign(obj.aux_k,1);
            assign(obj.aux_K,1);
            obj.X=X;
            obj.INIT_Y=X;
            obj.tris=tris;
            obj.C=C;
            obj.objType=objType;
            obj.spaceType=spaceType;
            obj.faces(size(tris,1))=Face();
            obj.source_dim=size(X,2);
            obj.target_dim=size(X,2);
            if exist('Y','var')
                obj.Y=Y;
            else
                obj.Y=sdpvar(size(X,1),obj.target_dim);
            end
            k = sdpvar(size(tris,1),1);
            K = sdpvar(size(tris,1),1);
            % calculate the differentials
            for i = 1:size(tris,1),
                tri=tris(i,:);
                s_tri=X(tri,:)';
                B = eye(size(s_tri,2)) - 1/size(s_tri,2);
                P = B/(s_tri*B);
                P(abs(P)<=1e2*eps) = 0;
                obj.faces(i).A = obj.Y(tri,:)'*P;
                obj.faces(i).K=K(i);
                obj.faces(i).k=k(i);
            end
            obj.log.t_Problem=stopper.stop;
            obj.setYalmipParams();
        end
        
        function setYalmipParams(obj)
            % set yalmip optimization parameters
            obj.yalmip_prm = sdpsettings;
            obj.yalmip_prm.solver = 'mosek';
            obj.yalmip_prm.verbose = 0;
        end
        
        function c=generateSpaceConstraints(obj)
            % generate constraints on singular values
            stopper=Stopper('Generating space constraints (%s)... ', char(obj.spaceType));
            c=cell(length(obj.faces)+1,1);
            if obj.spaceType==SpaceEnum.NONE
                return; % case of no constraints
            end
            % generate per-face constraints
            for i=1:length(obj.faces)
                face=obj.faces(i);
                A=face.A;
                F=face.F;
                K=face.K;
                k=face.k;
                if obj.objType==ObjectiveEnum.FEASIBILITY
                    K=K+obj.tau;
                    k=k-obj.tau;
                end
                FA = F'*A;
                constraint_upper = ([K*eye(size(A)) A; A' K*eye(size(A))]>=0); % bound maximal singular value (eq. (5) in the paper)
                constraint_lower = (((FA+FA')/2 - k*eye(size(A,1)))>=0); % bound the minimal singular value (eq. (6) and (11.d) in the paper)
                c{i}=constraint_upper+constraint_lower;
            end
            % add constraints over k,K (gamma, Gamma in the paper)
            if ~isa(obj.faces(1).K,'double') % do not try to constrain if K are doubles
                switch obj.spaceType
                    case SpaceEnum.BD % bounded conformal distortion
                        c{end}=([obj.faces.K]==obj.C*[obj.faces.k]);
                    case SpaceEnum.BI % bounded isometric distortion
                        c{end}=(([obj.faces.k]>=1/obj.C)+([obj.faces.K]<=obj.C));
                    case SpaceEnum.BSI % bounded scaled isometric distortion
                        c{end}=((obj.aux_k<=[obj.faces.k])+...
                            ([obj.faces.K]<=obj.aux_K)+...
                            (obj.aux_K==(obj.C^2)*obj.aux_k));
                    otherwise
                        error('constraint type not properly defined in code');
                end
            else
                warning('k,K seem to be constants -> k,K constraints are not generated');
            end
            % wrap up
            fprintf(' %d face constraints generated... ',nnz(~cellfun(@isempty,c)));
            obj.log.t_generateObjective=stopper.stop;
        end
        
        function [o, o_c]=generateObjective(obj)
            % generate objective for a mesh optimization problem
            stopper=Stopper('Generating objective (%s)... ', char(obj.objType));
            if obj.recalc_objective
                obj.objective = [];
                obj.objectiveConstraints = [];
                switch obj.objType % see ObjectiveEnum for more details
                    case ObjectiveEnum.ARAP
                        A = reshape([obj.faces.A],obj.target_dim,obj.target_dim,[]);
                        F = reshape([obj.faces.F],obj.target_dim,obj.target_dim,[]);
                        w = sqrt(obj.vols/sum(obj.vols));
                        W = bsxfun(@times,ones(size(A(:,:,1))),permute(w,[3 2 1]));
                        z = W(:).*(A(:)-F(:));
                        t = sdpvar;
                        obj.objective = t;
                        obj.objectiveConstraints = cone(z,t);
                    case ObjectiveEnum.FEASIBILITY
                        assert(obj.spaceType~=SpaceEnum.NONE);
                        obj.objective=obj.tau;
                    case ObjectiveEnum.DIRICHLET
                        A = reshape([obj.faces.A],obj.target_dim,obj.target_dim,[]);
                        w = sqrt(obj.vols/sum(obj.vols));
                        W = bsxfun(@times,ones(size(A(:,:,1))),permute(w,[3 2 1]));
                        z = A.*W;
                        z = z(:);
                        t = sdpvar;
                        obj.objective = t;
                        obj.objectiveConstraints = cone(z,t);
                    case ObjectiveEnum.AAAP
                        if obj.target_dim==3
                            % get adjacent faces
                            TR = triangulation(obj.tris,obj.X);
                            N = TR.neighbors;
                            I = repmat((1:size(N,1))',[1 4]);
                            N = N(:);
                            I = I(:);
                            ff = isnan(N) | (I>N);
                            N(ff) = [];
                            I(ff) = [];
                            % calc functional
                            A = reshape([obj.faces.A],obj.target_dim,obj.target_dim,[]);
                            Adiff = A(:,:,I)-A(:,:,N);
                            vols_sum = obj.vols(N)+obj.vols(I);
                            w = sqrt(vols_sum/sum(vols_sum));
                            W = bsxfun(@times,ones(size(Adiff(:,:,1))),permute(w,[3 2 1]));
                            z = W(:).*Adiff(:);
                            t_smooth = sdpvar;
                            obj.objective = t_smooth;
                            obj.objectiveConstraints = cone(z,t_smooth);
                            % need only be calculated once
                            obj.recalc_objective = false;
                        else
                            error('invalid dimension')
                        end
                    otherwise
                        error('no such functional');
                end
            else
                fprintf('SKIPPING ');
            end
            o=obj.objective;
            o_c=obj.objectiveConstraints;
            obj.log.t_generateSpaceConstraints=stopper.stop;
        end
        
        function updateDistortion(obj)
            % update the distortions
            Aval=reshape(double(cat(2,obj.faces.A)),size(obj.faces(1).A,1),size(obj.faces(1).A,2),[]);
            switch obj.spaceType
                case {SpaceEnum.BD, SpaceEnum.NONE}
                    for i = 1:length(obj.faces),
                        A = Aval(:,:,i);
                        obj.distortion(i) = cond(A);
                        obj.flipped(i) = det(A)<0;
                    end
                case SpaceEnum.BI
                    for i = 1:length(obj.faces),
                        A = Aval(:,:,i);
                        s = svd(A);
                        obj.distortion(i) = max(s(1),1/s(end));
                        obj.flipped(i) = det(A)<0;
                    end
                case SpaceEnum.BSI
                    distCtr = sqrt(double(obj.aux_k)*double(obj.aux_K));
                    for i = 1:length(obj.faces),
                        A = Aval(:,:,i);
                        s = svd(A./distCtr);
                        obj.distortion(i) = max(s(1),1/s(end));
                        obj.flipped(i) = det(A)<0;
                    end
                otherwise
                    error('Distortion is not defined');
            end
        end
        
        function updateFrames(obj)
            % update the rotations (frames) -- corresponding to Lemma 1 in the paper
            Aval=reshape(double(cat(2,obj.faces.A)),size(obj.faces(1).A,1),size(obj.faces(1).A,2),[]);
            for i = 1:length(obj.faces),
                A = Aval(:,:,i);
                [U, ~, V] = closestRotation(A);
                obj.faces(i).F = U*V';
            end
        end
        
        function c=gatherConstraints(obj,o_c,s_c,a_c)
            % collect the problem's constraints
            stopper=Stopper('Gathering constraints... ');
            if isa(a_c,'function_handle')
                a_c=a_c();
                a_c=[a_c{:}];
            end
            c=catRecursive(s_c{:})+o_c+a_c;
            obj.log.t_gatherConstraints=stopper.stop;
        end;
        
        function runSolver(obj,c,o)
            % use Yalmip to solve
            stopper=Stopper('Calling solvesdp... ');
            res = solvesdp(c,o,obj.yalmip_prm);
            if(res.problem)
                warning(yalmiperror(res.problem));
            end
            obj.objVal=double(o);
            obj.log.solverStatus=res.problem;
            obj.log.t_runSolver=stopper.stop;
            obj.log.t_solver=res.solvertime;
            fprintf('Solver time: (%.2g secs)\n', res.solvertime);
        end;
        
        function optimize(obj)
            % optimize - a single iteration of Algorithm 2
            stopper=Stopper('Optimizing:\n');
            [o, o_c]=obj.generateObjective();
            s_c=obj.generateSpaceConstraints();
            c=obj.gatherConstraints(o_c,s_c,obj.auxConstraints);
            obj.runSolver(c,o);
            fprintf('Overall optimization time: ');
            obj.log.t_optimize=stopper.stop;
        end
        
        function visualize(obj)
            % visualize results
            tri=obj.tris;
            Y=double(obj.Y);
            switch size(Y,2),
                case 2,
                    trimesh(tri,Y(:,1),Y(:,2));
                case 3,
                    tetramesh(tri,Y,'facealpha',0.1);
                otherwise,
                    error('unsupported dimension');
            end;
            d=max(obj.distortion);
            f=nnz(obj.flipped);
            title(sprintf('Distortion: %g\nFlipped: %d\nObjective: %d',d,f,obj.objVal));
            axis equal;
        end
        
        function dispReport(obj)
            % print report
            d=max(obj.distortion);
            f=nnz(obj.flipped);
            fprintf('Distortion: %g\nFlipped: %d\nObjective: %d\n',d,f,obj.objVal);
        end
        
        function setFrames(obj,source)
            % reset frames (used for initialization)
            if isempty(source)
                Id = eye(obj.target_dim,obj.target_dim);
                for i=1:length(obj.faces)
                    obj.faces(i).F = Id;
                end
            elseif isa(source,'Problem')
                [obj.faces.F] = source.faces.F;
            elseif isequal(size(source),size(obj.Y)),
                obj.INIT_Y=source;
                assign(obj.Y, source);
                obj.updateFrames;
                obj.updateDistortion;
            else
                error('invalid source type')
            end;
        end
        
        function calcObjectiveVal(obj)
            % calculate objective (useful when objective is defined in epigraph form
            try
                stopper=Stopper('Initializing objective:\n');
                [o, o_c]=obj.generateObjective();
                c = o_c + (obj.Y==double(obj.Y));
                solvesdp(c, o, obj.yalmip_prm);
                obj.objVal=double(o);
                obj.log.t_calcObjectiveVal=stopper.stop;
            catch
                warning('Objective value could not be initialized');
            end
        end
    end
end
