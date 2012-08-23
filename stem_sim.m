%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_sim < handle
    %stem simulator class
    
    properties
        stem_model=[];      %stem_model object
        nan_rate=[0 0];     %missing data rate percentage
        nan_pattern_par=[]; %spatial missing data pattern parameter theta. Missing data pattern is simulated by considering an exponential correlation function in the form exp(-d/theta) with d the euclidean distance between two sites.
    end
    
    methods
        function obj = stem_sim(stem_model)
            % constructor
            if nargin>=1
                if strcmp(class(stem_model),'stem_model')
                    obj.stem_model=stem_model;
                else
                    error('The input argument must be of class stem_model');
                end
            end
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            % simulate STEM3 data
            % see properties for details
            disp('Simulation started...');
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            if nargin>=2
                obj.nan_rate=nan_rate;
            else
                obj.nan_rate=[];
            end
            if nargin>=3
                obj.nan_pattern_par=nan_pattern_par;
            else
                obj.nan_pattern_par=[];
            end
            
            if not(isempty(obj.nan_rate))
                soglia_nan_g=norminv(1-obj.nan_rate(1)/2,0,1);
                nanmat_g=exp(-obj.stem_model.stem_data.DistMat_g./obj.nan_pattern_par(1));
            else
                nanmat_g=[];
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))&&(not(isempty(obj.nan_rate)))
                soglia_nan_r=norminv(1-obj.nan_rate(2)/2,0,1);
                nanmat_r=exp(-obj.stem_model.stem_data.DistMat_r./obj.nan_pattern_par(2));
            else
                nanmat_r=[];
            end
            nancov_g=[];
            if not(isempty(nanmat_g))
                for i=1:obj.stem_model.stem_data.stem_varset_g.nvar
                    nancov_g=blkdiag(nancov_g,get_block(obj.stem_model.stem_data.stem_varset_g.dim,i,obj.stem_model.stem_data.stem_varset_g.dim,i,nanmat_g));
                end
            end
            nancov_r=[];
            if not(isempty(nanmat_r))
                for i=1:obj.stem_model.stem_data.stem_varset_r.nvar
                    nancov_r=blkdiag(nancov_r,get_block(obj.stem_model.stem_data.stem_varset_r.dim,i,obj.stem_model.stem_data.stem_varset_r.dim,i,nanmat_r));
                end
            end            
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,j_rg,j_g] = obj.stem_model.get_sigma();

            if obj.stem_model.stem_par.p>0
                mu0=zeros(obj.stem_model.stem_par.p,1);
                sigma_0=0.01*eye(obj.stem_model.stem_par.p);
                Z=stem_sim.ar1_sim(obj.stem_model.stem_par.G,obj.stem_model.stem_par.sigma_eta,T,mu0,sigma_0);
            end
            
            Y=zeros(N,T);
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                if not(isempty(obj.stem_model.stem_data.stem_varset_r.X_rg))
                    W_r=mvnrnd(zeros(obj.stem_model.stem_data.stem_varset_r.N,1),sigma_W_r,T)';
                end
            end
            if obj.stem_model.stem_par.k>0
                W_g=zeros(obj.stem_model.stem_data.stem_varset_g.N,T,obj.stem_model.stem_par.k);
                for k=1:obj.stem_model.stem_par.k
                    W_g(:,:,k)=mvnrnd(zeros(obj.stem_model.stem_data.stem_varset_g.N,1),sigma_W_g{k},T)';
                end
            end
            W_eps=mvnrnd(zeros(N,1),sigma_eps,T)';
            
            for t=1:T
                if not(isempty(obj.stem_model.stem_data.X_rg))
                    if obj.stem_model.stem_data.X_rg_tv
                        Y(:,t)=Y(:,t)+D_apply(D_apply(M_apply(W_r(:,t),obj.stem_model.stem_data.M,'l'),obj.stem_model.stem_data.X_rg(:,1,t),'l'),j_rg,'l');
                    else
                        Y(:,t)=Y(:,t)+D_apply(D_apply(M_apply(W_r(:,t),obj.stem_model.stem_data.M,'l'),obj.stem_model.stem_data.X_rg(:,1,1),'l'),j_rg,'l');
                    end
                end
                if not(isempty(obj.stem_model.stem_data.X_beta))
                    if obj.stem_model.stem_data.X_beta_tv
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_beta(:,:,t)*obj.stem_model.stem_par.beta;
                    else
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_beta(:,:,1)*obj.stem_model.stem_par.beta;
                    end
                end
                if not(isempty(obj.stem_model.stem_data.X_time))
                    if obj.stem_model.stem_data.X_time_tv
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_time(:,:,t)*Z(:,t);
                    else
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_time(:,:,1)*Z(:,t);
                    end
                end      
                if not(isempty(obj.stem_model.stem_data.X_g))
                    if obj.stem_model.stem_data.X_g_tv
                        for k=1:obj.stem_model.stem_par.k
                            Y(:,t)=Y(:,t)+D_apply(D_apply(W_g(:,t,k),obj.stem_model.stem_data.X_g(:,1,t,k),'l'),j_g(:,k),'l');
                        end
                    else
                        for k=1:obj.stem_model.stem_par.k
                            Y(:,t)=Y(:,t)+D_apply(D_apply(W_g(:,t,k),obj.stem_model.stem_data.X_g(:,1,1,k),'l'),j_g(:,k),'l');    
                        end
                    end
                end
                Y(:,t)=Y(:,t)+W_eps(:,t);

                if not(isempty(nancov_g)&&isempty(nancov_r))
                    nanfill_g=mvnrnd(zeros(size(nancov_g,1),1),nancov_g);
                    nanfill_g(abs(nanfill_g)>=soglia_nan_g)=NaN;
                    nanfill_g(abs(nanfill_g)<soglia_nan_g)=1;
                    if not(isempty(nancov_r))
                        nanfill_r=mvnrnd(zeros(size(nancov_r,1),1),nancov_r);
                        nanfill_r(abs(nanfill_r)>=soglia_nan_r)=NaN;
                        nanfill_r(abs(nanfill_r)<soglia_nan_r)=1;
                    else
                        nanfill_r=[];
                    end
                    nanfill=[nanfill_g,nanfill_r];
                    Y(:,t)=Y(:,t).*nanfill';
                end
            end

            blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_g.dim)];
            Y_temp=[];
            for i=1:obj.stem_model.stem_data.stem_varset_g.nvar
                Y_temp{i}=Y(blocks(i)+1:blocks(i+1),:);
            end
            obj.stem_model.stem_data.stem_varset_g.Y=Y_temp;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                temp=max(blocks);
                blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_r.dim)]+temp;
                Y_temp=[];
                for i=1:obj.stem_model.stem_data.stem_varset_r.nvar
                    Y_temp{i}=Y(blocks(i)+1:blocks(i+1),:);
                end
                obj.stem_model.stem_data.stem_varset_r.Y=Y_temp;                
            end
            obj.stem_model.stem_data.update_data();
            obj.stem_model.stem_data.simulated=1;
            disp('Simulation ended. New Y updated.');
            disp('');
        end
        
        function set.nan_rate(obj,nan_rate)
            if not(isempty(nan_rate))
                if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                    if not(length(nan_rate)==2)
                        error('nan_rate must be a 2x1 vector');
                    end
                    for i=1:2
                        if (nan_rate(i)<0)||(nan_rate(i)>=1)
                            error('The nan_rate elements must be in the interval [0,1)');
                        end
                    end
                else
                    if not(isscalar(nan_rate))
                        error('nan_rate must be a scalar');
                    end
                    if (nan_rate<0)||(nan_rate>=1)
                        error('The nan_rate must be in the interval [0,1)');
                    end
                end
            end
            obj.nan_rate=nan_rate;
        end
        
        function set.nan_pattern_par(obj,nan_pattern_par)
            if not(isempty(nan_pattern_par))
                if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                    if not(length(nan_pattern_par)==2)
                        error('nan_pattern_par must be a 2x1 vector');
                    end
                    for i=1:2
                        if nan_pattern_par(i)<=0
                            error('The nan_pattern_par elements must be > 0');
                        end
                    end
                else
                    if not(isscalar(nan_pattern_par))
                        error('nan_rate must be a scalar');
                    end
                    if nan_pattern_par<=0
                        error('The nan_pattern_par must be between > 0');
                    end
                end
            end
            obj.nan_pattern_par=nan_pattern_par;
        end
    end
    
    methods (Static)
        function Z = ar1_sim(G,sigma_eta,T,mu0,sigma_0)
            Z=zeros(length(G),T+1);
            Z0=mvnrnd(mu0,sigma_0)';
            Z(:,1)=Z0;
            for t=2:T+1
                Z(:,t)=G*Z(:,t-1)+mvnrnd(zeros(length(G),1),sigma_eta)';
            end
            Z=Z(:,2:T+1);
        end
    end
    
end

