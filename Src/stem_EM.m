%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.


classdef stem_EM < EM
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    
    %CONSTANTS
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %S   = 2 if both point and pixel data are considered. S = 1 if only point data are considered.
    %T   - number of temporal steps
    %TT  = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    %p   - dimension of the latent temporal variable z
    
    properties
        stem_model=[];               %[stem_model object] (1x1) an object of class stem_model
        stem_EM_options=[];          %[stem_EM_options object] (1x1) an object of class stem_EM_options
    end
    
    methods
        function obj = stem_EM(stem_model,stem_EM_options)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_model        - [stem_model object] (1x1) an object of class stem_model
            %stem_EM_options   - [stem_EM_options]   (1x1) an object of class stem_EM_options
            %
            %OUTPUT
            %obj               - [stem_EM object]    (1x1) 
            if nargin<2
                error('All the arguments must be provided');
            end
            
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('The first argument must be of class stem_model');
            end
            
            if isa(stem_EM_options,'stem_EM_options')
                obj.stem_EM_options=stem_EM_options;
            else
                error('The second argument must be of class stem_EM_options');
            end
            
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
        end
        
        function st_EM_result = estimate(obj)
            %DESCRIPTION: EM estimation
            %
            %INPUT
            %obj            - [stem_EM object]      (1x1) the stem_EM object
            %
            %OUTPUT
            %st_EM_result   - [stem_EM_result object] (1x1) the stem_EM_result object
            
            t1_full=clock;
            if isempty(obj.stem_model)&&(nargin==0)
                error('You have to set the stem_model property first');
            end
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
            delta=9999;
            delta_logL=9999;
            last_logL=0;
            last_stem_par=obj.stem_model.par_vec;
            iteration=0;
            st_EM_result=stem_EM_result();
            st_EM_result.max_iterations=obj.stem_EM_options.max_iterations;
            st_EM_result.exit_toll=obj.stem_EM_options.exit_toll;
            st_EM_result.block_tapering_block_size=obj.stem_EM_options.block_tapering_block_size;
            st_EM_result.machine=computer;
            st_EM_result.date_start=datestr(now);
            exit_file=0;
            model_changed=0;
            while (delta>obj.stem_EM_options.exit_toll)&&(delta_logL>obj.stem_EM_options.exit_toll)&&(iteration<obj.stem_EM_options.max_iterations)&&(exit_file==0)||(model_changed==1)
                ct1=clock;
                iteration=iteration+1;
                disp('************************');
                disp(['Iteration ',num2str(iteration),' started...']);
                disp('************************');
                
                clear E_wb_y1
                clear sum_Var_wb_y1
                clear diag_Var_wb_y1
                clear cov_wb_z_y1
                clear E_wg_y1
                clear sum_Var_wp_y1
                clear diag_Var_wp_y1
                clear cov_wp_z_y1
                clear M_cov_wb_wp_y1
                clear cov_wpk_wph_y1
                clear diag_Var_e_y1
                clear E_e_y1
                clear sigma_eps
                clear Xbeta
                
                if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                   if not(isempty(obj.stem_model.stem_data.X_beta))
                       obj.stem_model.stem_data.X_beta=[];
                       nc=size(obj.stem_model.stem_data.stem_varset_p.X_beta{1},2)/2;
                       for t=1:size(obj.stem_model.stem_data.stem_varset_p.X_beta{1},3)
                           obj.stem_model.stem_data.X_beta(:,1:nc,t)=obj.stem_model.stem_data.stem_varset_p.X_beta{1}(:,1:nc,t).*obj.stem_model.stem_data.X_z;
                           obj.stem_model.stem_data.X_beta(:,nc+1:nc*2,t)=obj.stem_model.stem_data.stem_varset_p.X_beta{1}(:,nc+1:nc*2,t).*obj.stem_model.stem_data.X_z;
                       end
                   end
                end
                [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result] = obj.E_step();
                model_changed = obj.M_step(E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,iteration);

                if not(isempty(st_kalmansmoother_result))
                    if not(st_kalmansmoother_result.logL==0)
                        logL=st_kalmansmoother_result.logL;
                        st_EM_result.logL_all(iteration)=logL;
                        delta_logL=abs(logL-last_logL)/abs(logL);
                        last_logL=logL;
                        disp('****************');
                        disp( ['logL: ',num2str(logL)]);
                        disp(['logL relative delta: ',num2str(delta_logL)]);
                    else
                        delta_logL=9999;
                    end
                else
                    delta_logL=9999;
                end

                if (model_changed==0)
                    delta=norm(obj.stem_model.par_vec()-last_stem_par)/norm(last_stem_par);
                else
                    delta=9999;
                end
                
                last_stem_par=obj.stem_model.par_vec;
                disp(['Parameter delta norm: ',num2str(delta)]);
                obj.stem_model.print_par;
                ct2=clock;
                disp('**********************************************');
                disp(['Iteration ',num2str(iteration),' ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                disp('**********************************************');
                files=dir('exit.txt');
                if not(isempty(files))
                    exit_file=1;
                    delete('exit.txt');
                    disp('EM terminated because of the presence of the exit.txt file in the folder. The file exit.txt has been deleted.');
                end
            end
            

            npars = length(obj.stem_model.par_vec);
            if not(isempty(st_EM_result.logL_all))
                st_EM_result.AIC = st_EM_result.logL_all(end)-2*npars;
            end
                
            t2_full=clock;
            st_EM_result.stem_par=obj.stem_model.stem_par;
            st_EM_result.stem_kalmansmoother_result=st_kalmansmoother_result;
            st_EM_result.E_wp_y1=E_wp_y1;
            st_EM_result.E_wb_y1=E_wb_y1;
            st_EM_result.diag_Var_wp_y1=diag_Var_wp_y1;
            st_EM_result.diag_Var_wb_y1=diag_Var_wb_y1;
            
            y_hat=obj.stem_model.stem_data.Y;
            y_hat(isnan(y_hat))=0;
            y_hat=y_hat-E_e_y1;
            res=obj.stem_model.stem_data.Y-y_hat;
            
            mse=(nanvar(res'))'; % res is N*T matrix, as well as the obj.stem_model.stem_data.Y. Thus the R2 is for each site.
            st_EM_result.R2=1-(mse'./nanvar(obj.stem_model.stem_data.Y'))';
            
            blocks=[0 cumsum(obj.stem_model.dim)];
            counter=1;
            for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                st_EM_result.y_hat{counter}=y_hat(blocks(counter)+1:blocks(counter+1),:);
                st_EM_result.diag_Var_e_y1{counter}=diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                st_EM_result.res{counter}=res(blocks(counter)+1:blocks(counter+1),:);
                counter=counter+1;
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    st_EM_result.y_hat{counter}=y_hat(blocks(counter)+1:blocks(counter+1),:);
                    st_EM_result.diag_Var_e_y1{counter}=diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                    st_EM_result.res{counter}=res(blocks(counter)+1:blocks(counter+1),:);
                    counter=counter+1;
                end
            end
         
            blocks=[0 cumsum(obj.stem_model.dim)];
            counter=1;
            for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                y_hat_back=st_EM_result.y_hat{counter};
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                if obj.stem_model.stem_data.stem_varset_p.standardized
                    s=obj.stem_model.stem_data.stem_varset_p.Y_stds{i};
                    m=obj.stem_model.stem_data.stem_varset_p.Y_means{i};
                else
                    s=1;
                    m=0;
                end
                if not(obj.stem_model.stem_data.stem_varset_p.log_transformed)&&not(obj.stem_model.stem_data.stem_varset_p.boxcox_transformed)
                    y_hat_back=st_EM_result.y_hat{counter}*s+m;
                    y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                    var_y_hat=st_EM_result.diag_Var_e_y1{counter};
                    var_y_hat_back=var_y_hat*s^2;
                end
                if obj.stem_model.stem_data.stem_varset_p.log_transformed
                    y_hat_back=st_EM_result.y_hat{counter};
                    var_y_hat=st_EM_result.diag_Var_e_y1{counter};
                    y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                    y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                    var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat{counter}*s+m)+(var_y_hat*s^2));
                end
                
                if obj.stem_model.stem_data.stem_varset_p.boxcox_transformed
                    y_hat_back=st_EM_result.y_hat{counter};
                    var_y_hat=st_EM_result.diag_Var_e_y1{counter};
                    lambda=obj.stem_model.stem_data.stem_varset_p.Y_lambda{i};
                    
                    y_hat_back=((lambda*(y_hat_back*s+m)+1).^(1/lambda)).*(1+(var_y_hat*s^2*(1-lambda))./(2*(lambda*(y_hat_back*s+m)+1).^2));
                    y_back=(lambda*(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m)+1).^(1/lambda);
                    var_y_hat_back=0; %TO DO!
                end
                
                st_EM_result.y_hat_back{counter}=y_hat_back;
                st_EM_result.y_back{counter}=y_back;
                st_EM_result.res_back{counter}=y_back-y_hat_back;
                
                if not(isempty(var_y_hat_back))
                    st_EM_result.diag_Var_y_hat_back{counter} = var_y_hat_back;
                end
                counter=counter+1;
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                y_hat_back=st_EM_result.y_hat{counter};
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    if obj.stem_model.stem_data.stem_varset_b.standardized
                        s=obj.stem_model.stem_data.stem_varset_b.Y_stds{i};
                        m=obj.stem_model.stem_data.stem_varset_b.Y_means{i};
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&not(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        y_hat_back=st_EM_result.y_hat{counter}*s+m;
                        y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                        var_y_hat=st_EM_result.diag_Var_e_y1{counter};
                        var_y_hat_back=var_y_hat*s^2;
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        var_y_hat=st_EM_result.diag_Var_e_y1{counter};
                        y_hat_back=st_EM_result.y_hat{counter};
                        y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                        y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                        var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat{counter}*s+m)+(var_y_hat*s^2));
                    end
                    st_EM_result.y_hat_back{counter}=y_hat_back;
                    st_EM_result.y_back{counter}=y_back;
                   
                    if not(isempty(var_y_hat_back))
                        st_EM_result.diag_Var_y_hat_back{counter} = var_y_hat_back;
                    end
                    st_EM_result.res_back{counter}=y_back-y_hat_back;
                    counter=counter+1;
                end
            end
            
            st_EM_result.iterations=iteration;
            st_EM_result.computation_time=etime(t2_full,t1_full);
        end
        
        function [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result] = E_step(obj,T)
            %DESCRIPTION: E-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1) the stem_EM object
            %T                              - [integer >0]      (1x1) The E-step is computed only for the data related to the time steps between 1 and T
            %
            %OUTPUT
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp_k|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            %sigma_eps                      - [double]          (NxN) sigma_eps
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)
            
            if nargin==1
                T=obj.stem_model.stem_data.T;
            end
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            par=obj.stem_model.stem_par;
            
            disp('  E step started...');
            ct1_estep=clock;
          
            if p>0
                %Kalman smoother
                st_kalman=stem_kalman(obj.stem_model);
               
                [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,M] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps,0,[],[],obj.stem_EM_options.block_tapering_block_size,obj.stem_EM_options.workers);
                
                rr=size(sigma_Z,1);
                
                if not(obj.stem_model.stem_data.X_z_tv)
                    if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                        temp=obj.stem_model.stem_data.X_z{1};
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        if N-size(obj.stem_model.stem_data.X_z{1},1)>0
                            X_z_orlated=[obj.stem_model.stem_data.X_z{1};zeros(N-size(obj.stem_model.stem_data.X_z{1},1),size(obj.stem_model.stem_data.X_z{1},2))];
                        else
                            X_z_orlated=obj.stem_model.stem_data.X_z{1};
                        end
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                        if obj.stem_model.tapering
                            %is it possible to improve the sparse matrix var_Zt?
                            var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                        else
                            var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                        end
                    else
                        var_Zt=[];
                    end
                end
                if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                    var_Yt=sigma_geo+var_Zt;
                end
            else
                %[sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,~,~,aj_bp,aj_p,~,M] = obj.stem_model.get_sigma();
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,~,~,aj_bp,M] = obj.stem_model.get_sigma();

                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[],[],[]);
                var_Zt=[];
                rr=0;
                
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
            end
            
            aj_p=nan(Np+Nb,K);
            for k=1:K
                aj_p(:,k)=[ones(Np,1); zeros(Nb,1)];
            end
            
            if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                obj.stem_model.stem_data.stem_varset_p.Y{1}=[];
            end 
            
            E_e_y1=obj.stem_model.stem_data.Y;
            E_e_y1(isnan(E_e_y1))=0;
            if not(isempty(obj.stem_model.stem_data.X_beta))
                disp('    Xbeta evaluation started...');
                ct1=clock;
                Xbeta=zeros(N,T);
                if obj.stem_model.stem_data.X_beta_tv
                    for t=1:T
                        if size(obj.stem_model.stem_data.X_beta{t},1)<N
                            X_beta_orlated=[obj.stem_model.stem_data.X_beta{t};zeros(N-size(obj.stem_model.stem_data.X_beta{t},1),size(obj.stem_model.stem_data.X_beta{t},2))];
                        else
                            X_beta_orlated=obj.stem_model.stem_data.X_beta{t};
                        end
                        Xbeta(:,t)=X_beta_orlated*par.beta;
                    end
                else
                    if size(obj.stem_model.stem_data.X_beta{1},1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta{1};zeros(N-size(obj.stem_model.stem_data.X_beta{1},1),size(obj.stem_model.stem_data.X_beta{1},2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta{1};
                    end
                    Xbeta=repmat(X_beta_orlated*par.beta,1,T);
                end
                ct2=clock;
                disp(['    Xbeta evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                E_e_y1=E_e_y1-Xbeta;
            else
                Xbeta=[];
            end
            diag_Var_e_y1=zeros(N,T);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
            disp('    Conditional E, Var, Cov evaluation started...');
            ct1=clock;
            if not(isempty(obj.stem_model.stem_data.X_bp))
                if (obj.stem_model.tapering)
                    Lr=find(sigma_W_b);
                    nnz_b=length(Lr);
                end
                %cov_wb_yz time invariant case
                if not(obj.stem_model.stem_data.X_bp_tv)
                    cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp{1},'r'),aj_bp,'r');
                end
                E_wb_y1=zeros(Nb,T);
                if (obj.stem_model.tapering)
                    sum_Var_wb_y1=spalloc(size(sigma_W_b,1),size(sigma_W_b,2),nnz_b);
                else
                    sum_Var_wb_y1=zeros(Nb);
                end
                diag_Var_wb_y1=zeros(Nb,T);
                cov_wb_z_y1=zeros(Nb,rr,T);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_p))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_p{1});
                    nnz_p=length(Lg);
                end
                %cov_wp_yz time invariant case
                if not(obj.stem_model.stem_data.X_p_tv)
                    cov_wp_y=cell(K,1);
                    
                    for k=1:K
                        cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p{1}(:,k),'r'),aj_p(:,k),'r');
                    end
                end
                cov_wpk_wph_y1=cell(K,K);
                for h=1:K
                    for k=h+1:K
                        %these matrices can be sparse?
                        cov_wpk_wph_y1{k,h}=zeros(Np,T);
                    end
                end
                E_wp_y1=zeros(Np,T,K);
                sum_Var_wp_y1=cell(K,1);
                for k=1:K
                    if obj.stem_model.tapering
                        sum_Var_wp_y1{k}=spalloc(size(sigma_W_p{k},1),size(sigma_W_p{k},2),nnz_p);
                    else
                        sum_Var_wp_y1{k}=zeros(Np,Np);
                    end
                end
                diag_Var_wp_y1=zeros(Np,T,K);
                cov_wp_z_y1=zeros(Np,rr,T,K);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_bp)) && not(isempty(obj.stem_model.stem_data.X_p))
                M_cov_wb_wp_y1=zeros(N,T,K);
            else
                M_cov_wb_wp_y1=[];
            end
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                
                if obj.stem_model.stem_data.X_bp_tv
                    tBP=t;
                else
                    tBP=1;
                end
                if obj.stem_model.stem_data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if obj.stem_model.stem_data.X_p_tv
                    tP=t;
                else
                    tP=1;
                end
                
                %evaluate var_yt in the time variant case
                if obj.stem_model.stem_data.X_tv
                    if not(isempty(obj.stem_model.stem_data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_model.stem_data.X_bp{tBP},'b'),aj_bp,'b');
                    end
                    if not(isempty(obj.stem_model.stem_data.X_p))
                        if isempty(sigma_geo)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(obj.stem_model.stem_data.X_p{1},2)
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p{tP}(:,k),'b'),aj_p(:,k),'b');
                        end
                    end
                    if isempty(sigma_geo)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if p>0
                        if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                            temp=obj.stem_model.stem_data.X_z{tT};
                            temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[obj.stem_model.stem_data.X_z{tT};zeros(N-size(obj.stem_model.stem_data.X_z{tT},1),size(obj.stem_model.stem_data.X_z{tT},2))];
                        end

                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            if obj.stem_model.tapering
                                %is it possible to improve the sparse matrix var_Zt?
                                var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                            else
                                var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                            end
                            if isempty(sigma_geo)
                                var_Yt=var_Zt;
                            else
                                var_Yt=sigma_geo+var_Zt;
                            end
                        else
                            var_Zt=[];
                            var_Yt=[];
                        end
                    else
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            var_Yt=sigma_geo;
                        else
                            var_Yt=[];
                        end
                    end
                end
                
                if p>0
                    if obj.stem_model.product_step>0
                        blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                        end
                        if sum(obj.stem_EM_options.block_tapering_block_size)>0
                            Pk_s_sparse=sparse(st_kalmansmoother_result.Pk_s{t+1});
                        end
                        for i=1:length(blocks)-1
                            if sum(obj.stem_EM_options.block_tapering_block_size)>0
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(X_z_orlated(blocks(i)+1:blocks(i+1),:)*Pk_s_sparse*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                            else
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(X_z_orlated(blocks(i)+1:blocks(i+1),:)*st_kalmansmoother_result.Pk_s{t+1}*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                            end
                        end
                    else
                        temp=X_z_orlated*st_kalmansmoother_result.Pk_s{t+1};
                        diag_Var_e_y1(:,t)=diag(temp*X_z_orlated');
                    end
                    %update E(e|y1)
                    temp=st_kalmansmoother_result.zk_s(:,t+1);
                    E_e_y1(:,t)=E_e_y1(:,t)-X_z_orlated*temp;
                end
                
                if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), X_z_orlated(Lt,:)*sigma_Z; sigma_Z*X_z_orlated(Lt,:)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    
                    res=obj.stem_model.stem_data.Y;
                    if not(isempty(Xbeta))
                        res=res-Xbeta;
                    end
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=stem_misc.chol_solve(chol_H1t,temp2(r));
                    else
                        chol_H1t=chol(H1t);
                        cs=stem_misc.chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end

                if not(isempty(obj.stem_model.stem_data.X_bp))
                    %cov_wb_yz time variant case
                    if obj.stem_model.stem_data.X_bp_tv
                        cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp{tBP},'r'),aj_bp,'r');
                    end
                    cov_wb_y1z=[cov_wb_y(:,Lt),zeros(size(cov_wb_y,1),rr)];
                    %compute E(w_b|y1);
                    E_wb_y1(:,t)=cov_wb_y1z*cs;
                    %compute Var(w_b|y1)
                    if obj.stem_model.tapering
                        temp_b(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wb_y1z(:,r)');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    else
                        temp_b=stem_misc.chol_solve(chol_H1t,cov_wb_y1z');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    end
                    
                    if p>0
                        %compute cov(w_b,z|y1)
                        cov_wb_z_y1(:,:,t)=temp_b(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s{t+1};
                        Var_wb_y1=Var_wb_y1+cov_wb_z_y1(:,:,t)*temp_b(end-rr+1:end,:);
                        %update diag(Var(e|y1))
                        temp=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp{tBP},'l'),aj_bp,'l');
                        
                        if obj.stem_model.product_step>0
                            blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                        end
                    else
                        cov_wb_z_y1=[];
                    end
                    %compute diag(Var(w_b|y1))
                    diag_Var_wb_y1(:,t)=diag(Var_wb_y1);
                    %compute sum(Var(w_b|y1))
                    sum_Var_wb_y1=sum_Var_wb_y1+Var_wb_y1;
                    %update E(e|y1)
                    E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp{tBP},'l'),aj_bp,'l');
                    
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(diag_Var_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp{tBP},'b'),aj_bp,'b');
               
                else
                    E_wb_y1=[];
                    diag_Var_wb_y1=[];
                    sum_Var_wb_y1=[];
                    cov_wb_z_y1=[];
                end
                clear temp_b
                
                if not(isempty(obj.stem_model.stem_data.X_p))
                    if obj.stem_model.stem_data.X_p_tv
                        %cov_wp_yz time invariant case
                        cov_wp_y=cell(K,1);
                        for k=1:K
                            cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p{tP}(:,k),'r'),aj_p(:,k),'r');
                        end
                    end
                    temp_p=cell(K,1);
                    for k=1:K
                        cov_wp_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                        %compute E(w_p_k|y1);
                        E_wp_y1(:,t,k)=cov_wp_y1z*cs;
                        %compute Var(w_p_k|y1)
                        if obj.stem_model.tapering
                            temp_p{k}(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wp_y1z(:,r)');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        else
                            temp_p{k}=stem_misc.chol_solve(chol_H1t,cov_wp_y1z');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        end
                        
                        if p>0
                            %compute cov(w_p,z|y1)
                            cov_wp_z_y1(:,:,t,k)=temp_p{k}(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s{t+1};
                            Var_wp_y1=Var_wp_y1+cov_wp_z_y1(:,:,t,k)*temp_p{k}(end-rr+1:end,:);
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),obj.stem_model.stem_data.X_p{tP}(:,k),'l'),aj_p(:,k),'l');
                            if obj.stem_model.product_step>0
                                blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                            end
                        else
                            cov_wp_z_y1=[];
                        end
                        diag_Var_wp_y1(:,t,k)=diag(Var_wp_y1);
                        sum_Var_wp_y1{k}=sum_Var_wp_y1{k}+Var_wp_y1;
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,k),obj.stem_model.stem_data.X_p{tP}(:,k),'l'),aj_p(:,k),'l');

                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(diag_Var_wp_y1(:,t,k),obj.stem_model.stem_data.X_p{tP}(:,k),'b'),aj_p(:,k),'b'); %K varianze
                        
                        if not(isempty(obj.stem_model.stem_data.X_bp))
                            %compute M_cov(w_b,w_p|y1) namely M*cov(w_b,w_p|y1)
                            if obj.stem_model.product_step>0
                                blocks=0:obj.stem_model.product_step:length(M);
                                if not(blocks(end)==length(M))
                                    blocks=cat(2,blocks,length(M));
                                end
                                for i=1:length(blocks)-1
                                    if p>0
                                        M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1))+cov_wb_z_y1(M(blocks(i)+1:blocks(i+1)),:,t)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1))); %ha gia' l'stem_misc.M_apply su left!!
                                    else
                                        M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M))+cov_wb_z_y1(M,:,t)*temp_p{k}(end-rr+1:end,1:length(M))); %ha giï¿??? l'stem_misc.M_apply su left!!
                                else
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp{tBP},'l'),aj_bp,'l');
                            temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p{tP}(:,k);zeros(Nb,1)],'l'),aj_p(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_pk,w_ph|y1);
                        cov_wpk_wph_y1=cell(K,K);
                        for h=1:K
                            for k=h+1:K
                                cov_wpk_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                                if obj.stem_model.product_step>0
                                    blocks=0:obj.stem_model.product_step:size(cov_wpk_y1z,1);
                                    if not(blocks(end)==size(cov_wpk_y1z,1))
                                        blocks=cat(2,blocks,size(cov_wpk_y1z,1));
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wp_z_y1))
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1))+cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{h}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wp_z_y1))
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h}+cov_wp_z_y1(:,:,t,k)*temp_p{h}(end-rr+1:end,:));
                                    else
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h});
                                    end
                                end
                                temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),obj.stem_model.stem_data.X_p{tP}(:,k),'l'),aj_p(:,k),'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p{tP}(:,h);zeros(Nb,1)],'l'),aj_p(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    else
                        cov_wpk_wph_y1=[];
                    end
                else
                    E_wp_y1=[];
                    diag_Var_wp_y1=[];
                    sum_Var_wp_y1=[];
                    M_cov_wb_wp_y1=[];
                    cov_wp_z_y1=[];
                    cov_wpk_wph_y1=[];
                end
                %delete the variables the dimension of which changes every t
                clear temp_p
                clear temp
                if obj.stem_model.stem_data.X_tv
                    sigma_geo=[];
                end
            end
            
            if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                obj.stem_model.stem_data.stem_varset_p.Y{1}=obj.stem_model.stem_data.Y;
            end
            
            ct2=clock;
            disp(['    Conditional E, Var, Cov evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            ct2_estep=clock;
            disp(['  E step ended in ',stem_misc.decode_time(etime(ct2_estep,ct1_estep))]);
            disp('');
        end
        
        function model_changed = M_step(obj,E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,iteration)
            %DESCRIPTION: M-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1) the stem_EM object
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            %sigma_eps                      - [double]          (NxN|NxNxT) sigma_eps
            %st_kalmansmoother_result       - [st_kalmansmoother_result object] (1x1)
            %iteration                      - [double]          (1x1) EM iteration number
            %model_changed                  - [boolean]         (1x1) 1: the number of model parameters changed from the previous iteration (can happen with clustering); 0: no change
            %
            %OUTPUT
            %none: the stem_par property of the stem_model object is updated
            
            disp('  M step started...');
            ct1_mstep=clock;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            N=obj.stem_model.stem_data.N;
            T=obj.stem_model.stem_data.T;
            K=obj.stem_model.stem_par.k;
            M=obj.stem_model.stem_data.M;
            dim=obj.stem_model.stem_data.dim;
            
            par=obj.stem_model.stem_par;
            st_par_em_step=par;
            fminsearch_max_iter=obj.stem_EM_options.fminsearch_max_iterations;
            
            aj_bp=obj.stem_model.get_aj();
            
            if not(iscell(sigma_eps)) 
                if size(sigma_eps,3)==1
                    d=1./diag(sigma_eps);
                    I=1:length(d);
                    inv_sigma_eps=sparse(I,I,d);
                end
            end
            
            if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                obj.stem_model.stem_data.stem_varset_p.Y{1}=[];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             beta update                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_beta))
                ct1=clock;
                disp('    beta updating started...');
                temp1=zeros(size(obj.stem_model.stem_data.X_beta{1},2));
                temp2=zeros(size(obj.stem_model.stem_data.X_beta{1},2),1);
                if not(iscell(sigma_eps)) 
                    if size(sigma_eps,3)==1
                        d=diag(inv_sigma_eps);
                    end
                end
                
                for t=1:T
                    if iscell(sigma_eps) 
                        d=1./diag(sigma_eps{t});
                        I=1:length(d);
                        inv_sigma_eps=sparse(I,I,d);
                        d=diag(inv_sigma_eps);
                    else
                        if size(sigma_eps,3)>1
                            d=1./diag(sigma_eps(:,:,t));
                            I=1:length(d);
                            inv_sigma_eps=sparse(I,I,d);
                            d=diag(inv_sigma_eps);
                        end
                    end
                    
                    Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                    if obj.stem_model.stem_data.X_beta_tv
                        tB=t;
                    else
                        tB=1;
                    end
                    if size(obj.stem_model.stem_data.X_beta{tB},1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta{tB};zeros(N-size(obj.stem_model.stem_data.X_beta{tB},1),size(obj.stem_model.stem_data.X_beta{tB},2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta{tB};
                    end
                    temp1=temp1+X_beta_orlated(Lt,:)'*stem_misc.D_apply(X_beta_orlated(Lt,:),d(Lt),'l');
                    temp2=temp2+X_beta_orlated(Lt,:)'*stem_misc.D_apply(E_e_y1(Lt,t)+X_beta_orlated(Lt,:)*par.beta,d(Lt),'l');
                end
                st_par_em_step.beta=temp1\temp2;
                ct2=clock;
                disp(['    beta updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              sigma_eps                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if par.stem_par_constraints.sigma_eps_fixed==0
                disp('    sigma_eps updating started...');
                ct1=clock;
                
                if iscell(sigma_eps) 
                    if par.flag_logsigma==1
                        temp = [];
                        temp1 =zeros(N,1);
                        bmat = zeros(N*T,obj.stem_model.stem_par.k_sigma);
                        for t=1:T
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                            temp=cat(1,temp,temp1);
                            bmat(N*(t-1)+(1:N),:) = full(getbasismatrix(obj.stem_model.stem_data.X_h(:,t),...
                            obj.stem_model.stem_data.stem_fda.spline_basis_sigma));
                        end 
                        %ct1 =clock;
                        if obj.stem_EM_options.workers==1
                            for i =1:length(st_par_em_step.sigma_eps)
                                st_par_em_step.sigma_eps(i) = fminsearch(@(x) stem_misc.update_coe_sigma_eps(x,i,st_par_em_step.sigma_eps,temp,bmat),par.sigma_eps(i),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.001,'UseParallel','always'));
                            end
                        else
                            parfor i =1:length(st_par_em_step.sigma_eps)
                                sigma_eps_coe(i) = fminsearch(@(x) stem_misc.update_coe_sigma_eps(x,i,par.sigma_eps,temp,bmat),par.sigma_eps(i),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.001,'UseParallel','always'));
                            end
                            for i =1:length(st_par_em_step.sigma_eps)
                                st_par_em_step.sigma_eps(i)=sigma_eps_coe(i);
                            end
                        end
                    else
                        error('Code below must be tested first');
                        temp=zeros(par.k_sigma,1);
                        temp1=zeros(N,1);
                        temp2=zeros(par.k_sigma, par.k_sigma);
                        for t=1:T
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                            if sum(Lt)<N
                                temp1(~Lt)=full(getbasismatrix(obj.stem_model.stem_data.X_h(~Lt,t),...
                                    obj.stem_model.stem_data.stem_fda.spline_basis_sigma))*par.sigma_eps;
                            end
                            bmat = full(getbasismatrix(obj.stem_model.stem_data.X_h(:,t),...
                                obj.stem_model.stem_data.stem_fda.spline_basis_sigma));
                            temp = temp + bmat'*temp1;
                            temp2 = temp2 + bmat'*bmat;
                        end
                        st_par_em_step.sigma_eps = temp2\temp;
                    end
                else
                    if size(sigma_eps,3)==1
                        temp=zeros(N,1);
                        temp1=zeros(N,1);
                        d=diag(sigma_eps);
                        for t=1:T
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            %the next two lines are ok only for sigma_eps diagonal
                            temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                            temp1(~Lt)=d(~Lt);
                            temp=temp+temp1;
                        end
                        temp=temp/T;
                    else
                        temp=zeros(N,T);
                        temp1=zeros(N,1);
                        
                        for t=1:T
                            d=diag(sigma_eps(:,:,t));
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            %the next two lines are ok only for sigma_eps diagonal
                            temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                            temp1(~Lt)=d(~Lt);
                            temp(:,t)=temp1;
                        end
                    end

                    if not(obj.stem_model.stem_data.stem_modeltype.is('MBC'))
                        blocks=[0 cumsum(dim)];
                        if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))
                            if size(sigma_eps,3)==1
                                for i=1:length(dim)
                                    st_par_em_step.sigma_eps(i,i)=mean(temp(blocks(i)+1:blocks(i+1)));
                                end
                            else
                                for t=1:size(sigma_eps,3)
                                    for i=1:length(dim)
                                        st_par_em_step.sigma_eps(i,i,t)=mean(temp(blocks(i)+1:blocks(i+1),t));
                                    end
                                end
                            end
                        else
                            st_par_em_step.sigma_eps=mean(temp);
                        end
                    else
                        if obj.stem_model.stem_data.stem_modeltype.clustering_error_type_is('Shared')
                            st_par_em_step.sigma_eps(1,1)=mean(temp);
                        end
                        if obj.stem_model.stem_data.stem_modeltype.clustering_error_type_is('Proportional')
                            std_variance=obj.stem_model.stem_data.stem_varset_p.Y_stds{1}.^2;
                            std_variance=1./std_variance;
                            std_variance=std_variance./max(std_variance);
                            st_par_em_step.sigma_eps(1,1)=mean(temp./std_variance);
                        end
                        if obj.stem_model.stem_data.stem_modeltype.clustering_error_type_is('Dynamic')
                            [~,idx]=max(abs(obj.stem_model.stem_data.X_z{1}),[],2);
                            for i=1:obj.stem_model.stem_par.p
                                L=idx==i;
                                if sum(L)>0
                                    st_par_em_step.sigma_eps(i,i)=mean(temp(L));
                                end
                            end
                        end
                    end
                end

                ct2=clock;
                disp(['    sigma_eps updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %    G and sigma_eta    %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if par.p>0
                if not(obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    disp('    G and sigma_eta updating started...');
                    ct1=clock;
                    
                    if obj.stem_model.stem_data.stem_datestamp.irregular==0
                        %regular time steps case
                        if not(obj.stem_model.stem_par.stem_par_constraints.time_diagonal)
                            S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(cat(3,st_kalmansmoother_result.Pk_s{2:end}),3);
                            S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(cat(3,st_kalmansmoother_result.Pk_s{1:end-1}),3);
                            S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(cat(3,st_kalmansmoother_result.PPk_s{2:end}),3);
                        else
                            S11=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'))+diag(diag(sum(cat(3,st_kalmansmoother_result.Pk_s{2:end}),3)));
                            S00=diag(diag(st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(cat(3,st_kalmansmoother_result.Pk_s{1:end-1}),3)));
                            S10=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(cat(3,st_kalmansmoother_result.PPk_s{2:end}),3)));
                        end
                        
                        temp=abs(par.G-eye(size(par.G)));
                        if sum(temp(:))==0
                            %random walk
                            st_par_em_step.G=par.G;
                        else
                            temp=S10/S00;
                            if max(eig(temp))<1
                                st_par_em_step.G=temp;
                            else
                                warning('G is not stable. The last G is retained.');
                            end
                        end
                        
                        temp=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                        if min(eig(temp))>0
                            st_par_em_step.sigma_eta=temp;
                        else
                            warning('Sigma eta is not s.d.p. The last s.d.p. solution is retained');
                        end
                    else
                        %G
                        G_temp=par.G;

                        %first the diagonal elements are estimated
                        for h=1:size(G_temp,1)
                            initial=G_temp(h,h);
                            ctv1=clock;
                            min_result = fminsearch(@(x) stem_EM.function_g_element(x,h,h,G_temp,par.sigma_eta,par.lambda,st_kalmansmoother_result),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            ctv2=clock;
                            disp(['    G(',num2str(h),',',num2str(h),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            G_temp(h,h)=min_result;
                        end
                        if not(obj.stem_model.stem_par.stem_par_constraints.time_diagonal)
                            %if not diagonal, also the other elements are estimated
                            for h=1:size(G_temp,1)
                                for k=1:size(G_temp,1)
                                    if not(h==k)
                                        initial=G_temp(h,k);
                                        ctv1=clock;
                                        min_result = fminsearch(@(x) stem_EM.function_g_element(x,h,k,G_temp,par.sigma_eta,par.lambda,st_kalmansmoother_result),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                        ctv2=clock;
                                        disp(['    G(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                                        G_temp(h,k)=min_result;
                                    end
                                end
                            end
                        end
                        st_par_em_step.G=G_temp;
                        %when irregular, what's the meaning of stamp? it
                        %should provide as integer.
                        %sigma_eta
                        times=st_kalmansmoother_result.stem_datestamp.stamp;
                        %set the time of the state at t_0
                        m=mean(diff(times));
                        times=[times(1)-m,times];
                        zk_s=st_kalmansmoother_result.zk_s;
                        Pk_s=st_kalmansmoother_result.Pk_s;
                        PPk_s=st_kalmansmoother_result.PPk_s;
                        
                        U=cell(length(zk_s),1);
                        for t=2:length(zk_s)
                            delta=times(t)-times(t-1);
                            n_ref=delta/par.lambda;
                            G_t=st_par_em_step.G^n_ref;
                            U{t}=(zk_s(:,t)*zk_s(:,t)'+Pk_s{t})-((zk_s(:,t)*zk_s(:,t-1)'+PPk_s{t})*G_t')-...
                                (G_t*(zk_s(:,t)*zk_s(:,t-1)'+PPk_s{t})')+(G_t*(zk_s(:,t-1)*zk_s(:,t-1)'+Pk_s{t-1})*G_t');
                        end
                        
                        sigma_eta_temp=par.sigma_eta;
                        %first the variances are estimated
                        for k=1:size(sigma_eta_temp,1)
                            initial=sigma_eta_temp(k,k);
                            ctv1=clock;
                            min_result = fminsearch(@(x) stem_EM.function_sigma_eta_element(x,k,k,st_par_em_step.G,sigma_eta_temp,par.lambda,times,U),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            ctv2=clock;
                            disp(['    sigma_eta(',num2str(k),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            sigma_eta_temp(k,k)=min_result;
                        end
                        if not(obj.stem_model.stem_par.stem_par_constraints.time_diagonal)
                            %then the extra diagonal elements
                            for h=1:size(sigma_eta_temp,1)
                                for k=h+1:size(sigma_eta_temp,2)
                                    initial=sigma_eta_temp(h,k);
                                    ctv1=clock;
                                    min_result = fminsearch(@(x) stem_EM.function_sigma_eta_element(x,h,k,st_par_em_step.G,sigma_eta_temp,par.lambda,times,U),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    ctv2=clock;
                                    disp(['    sigma_eta(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                                    sigma_eta_temp(h,k)=min_result;
                                    sigma_eta_temp(k,h)=min_result;
                                end
                            end
                        end

                        if min(eig(sigma_eta_temp))>=0
                            st_par_em_step.sigma_eta=sigma_eta_temp;
                        else
                            disp('    sigma_eta is not positive definited. The last matrix is retained.');
                        end
                    end
                    ct2=clock;
                    disp(['    G and sigma_eta updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                else
                    if T>1
                        disp('    G updating started...');
                        ct1=clock;
                        
                        s=st_kalmansmoother_result.Pk_s{2};
                        for t=3:length(st_kalmansmoother_result.Pk_s)
                            s=s+st_kalmansmoother_result.Pk_s{t};
                        end
                        S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+s;
                        
                        s=st_kalmansmoother_result.Pk_s{1};
                        for t=2:length(st_kalmansmoother_result.Pk_s)-1
                            s=s+st_kalmansmoother_result.Pk_s{t};
                        end
                        S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+s;
                        
                        s=st_kalmansmoother_result.PPk_s{2};
                        for t=3:length(st_kalmansmoother_result.PPk_s)
                            s=s+st_kalmansmoother_result.PPk_s{t};
                        end
                        S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+s;
                        
                        temp=abs(par.G-eye(size(par.G)));
                        if sum(temp(:))==0
                            st_par_em_step.G=par.G;
                        else
                            temp=zeros(size(par.G));
                            if (par.p==par.q)
                                for i=1:par.p
                                    temp(i,i)=trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S10))/trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S00));
                                end
                            else
                                for i=1:par.p
                                    temp(i,i)=trace(stem_misc.get_block(ones(1,par.p)*dim(1),i,ones(1,par.p)*dim(1),i,S10))/trace(stem_misc.get_block(ones(1,par.p)*dim(1),i,ones(1,par.p)*dim(1),i,S00));
                                end
                            end
                            if max(eig(temp))<1
                                st_par_em_step.G=temp;
                            else
                                warning('G is not stable. The last G is retained.');
                            end
                            ct2=clock;
                            disp(['    G updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                        end
                        
                        dim=obj.stem_model.stem_data.dim;
                        if obj.stem_model.stem_data.stem_modeltype.is('f-HDGM')
                            G_tilde_diag=kron(diag(par.G),ones(dim(1),1));
                        else
                            G_tilde_diag=[];
                            for i=1:par.p
                                G_tilde_diag=cat(1,G_tilde_diag,par.G(i,i)*ones(dim(i),1));
                            end
                        end
                        G_tilde=sparse(1:length(G_tilde_diag),1:length(G_tilde_diag),G_tilde_diag,length(G_tilde_diag),length(G_tilde_diag));
                        temp=S11-S10*G_tilde'-G_tilde*S10'+G_tilde*S00*G_tilde';
                    else
                        temp=st_kalmansmoother_result.zk_s(:,2)*st_kalmansmoother_result.zk_s(:,2)'+st_kalmansmoother_result.Pk_s{2};
                    end

                    if obj.stem_model.stem_data.stem_modeltype.is('f-HDGM')
                        disp('    v_z and theta_z updating started...');
                        ct1=clock;
                        v_temp=par.v_z;
 
                        %indices are permutated in order to avoid deadlock
                        initial=cell(par.p);
                        if obj.stem_EM_options.workers==1    
                            kindex=randperm(size(par.v_z,1));
                            for j=kindex
                                initial{j}=par.v_z(j,j);
                                if not(isempty(obj.stem_model.stem_data.stem_gridlist_p.tap))
                                    error('Block tapering cannot be used when tapering is already enabled');
                                else
                                    dim=obj.stem_model.stem_data.stem_varset_p.dim;
                                    dim=ones(1,par.p)*dim(1);
                                    if sum(obj.stem_EM_options.block_tapering_block_size)==0
                                        [U,~,~] = stem_misc.get_block(dim,j,dim,j,temp);
                                        if not(strcmp(par.correlation_type,'expsphere'))
                                            [B,~,~] = stem_misc.get_block(dim,j,dim,j,obj.stem_model.DistMat_z);
                                        else
                                            for b=1:2
                                                [B{b},~,~] = stem_misc.get_block(dim,j,dim,j,obj.stem_model.DistMat_z{b});
                                            end
                                        end
                                        sigma_theta=stem_misc.correlation_function(par.theta_z(:,j),B,par.correlation_type);
                                        c=chol(sigma_theta);
                                        min_result{j}=trace(stem_misc.chol_solve(c,U))/T/size(c,1);
                                    else
                                        block_size=obj.stem_EM_options.block_tapering_block_size;
                                        min_result{j} = fminsearch(@(x) stem_EM.geo_coreg_function_velement_block(x,j,j,par.v_z,par.theta_z(:,j),par.correlation_type,obj.stem_model.DistMat_z,...
                                            dim,temp,T,block_size),initial{j},optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                end
                                v_temp(j,j)=min_result{j};
                            end
                        else
                            tap=obj.stem_model.stem_data.stem_gridlist_p.tap;
                            p=par.p;
                            dim=ones(1,p)*obj.stem_model.stem_data.stem_varset_p.dim(1);
                            v_z=par.v_z;
                            correlation_type=par.correlation_type;
                            block_tapering_block_size=obj.stem_EM_options.block_tapering_block_size;
                            theta_z=par.theta_z;
                            DistMat_z=obj.stem_model.DistMat_z;
                            
                            parfor j=1:par.p
                                initial=v_z(j,j);
                                if not(isempty(tap))
                                    error('Block tapering cannot be used when tapering is already enabled');
                                else
                                    if sum(block_tapering_block_size)==0
                                        [U,~,~] = stem_misc.get_block(dim,j,dim,j,temp);
                                        if not(strcmp(correlation_type,'expsphere'))
                                            [B,~,~] = stem_misc.get_block(dim,j,dim,j,DistMat_z);
                                        else
                                            B=cell(2,1);
                                            for b=1:2
                                                [B{b},~,~] = stem_misc.get_block(dim,j,dim,j,DistMat_z{b});
                                            end
                                        end
                                        sigma_theta=stem_misc.correlation_function(theta_z(:,j),B,correlation_type);
                                        c=chol(sigma_theta);
                                        min_result{j}=trace(stem_misc.chol_solve(c,U))/T/size(c,1);
                                        
                                    else
                                        block_size=block_tapering_block_size;
                                        min_result{j} = fminsearch(@(x) stem_EM.geo_coreg_function_velement_block(x,j,j,v_z,theta_z(:,j),correlation_type,DistMat_z,...
                                            dim,temp,T,block_size),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                end
                            end
                            for j=1:par.p
                                v_temp(j,j)=min_result{j};
                            end
                        end

                        ct2=clock;
                        disp(['    v_z updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                        
                        if min(eig(v_temp))>=0
                            st_par_em_step.v_z=v_temp;
                        else
                            disp('    v_z is not positive definited. The last matrix is retained.');
                        end
                    end
                    
                                        
                    if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                        disp('    v_z and theta_z updating started...');
                        ct1=clock;
                        v_temp=par.v_z;
                        
                        %indices are permutated in order to avoid deadlock
                        if obj.stem_EM_options.workers==1
                            %indices are permutated in order to avoid deadlock
                            kindex=randperm(size(par.v_z,1));
                            for k=kindex
                                hindex=randperm(size(par.v_z,1)-k+1)+k-1;
                                for h=hindex
                                    initial=par.v_z(k,h);
                                    ctv1=clock;
                                    if sum(obj.stem_EM_options.block_tapering_block_size)==0
                                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    else
                                        if not(isempty(obj.stem_model.stem_data.stem_gridlist_p.tap))
                                            error('Block tapering cannot be used when tapering is already enabled');
                                        else
                                            block_size=obj.stem_EM_options.block_tapering_block_size;
                                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement_block(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,block_size),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                        end
                                    end
                                    ctv2=clock;
                                    disp(['    v_z(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                                    v_temp(k,h)=min_result;
                                    v_temp(h,k)=min_result;
                                end
                            end
                        else
                            Idx = [];
                            for h=1:par.p
                                for k=h:par.p
                                    Idx = cat(1,Idx,sub2ind(size(par.v_z),h,k));
                                end
                            end
                            parfor j=1:numel(Idx)
                                [h,k]=ind2sub(size(par.v_z),Idx(j));
                                initial{j}=par.v_z(k,h);
                                if sum(obj.stem_EM_options.block_tapering_block_size)==0
                                    min_result{j} = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                        obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial{j},optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                else
                                    if not(isempty(obj.stem_model.stem_data.stem_gridlist_p.tap))
                                        error('Block tapering cannot be used when tapering is already enabled');
                                    else
                                        block_size=obj.stem_EM_options.block_tapering_block_size;
                                        min_result{j} = fminsearch(@(x) stem_EM.geo_coreg_function_velement_block(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,block_size),initial{j},optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                end
                            end
                            for i=1:numel(Idx)
                                [h,k]=ind2sub(size(par.v_z),Idx(i));
                                v_temp(k,h)=min_result{i};
                                v_temp(h,k)=min_result{i};
                            end
                        end
                        ct2=clock;
                        disp(['    v_z and theta_z updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                        
                        if min(eig(v_temp))>=0
                            st_par_em_step.v_z=v_temp;
                        else
                            disp('    v_z is not positive definited. The last matrix is retained.');
                        end
                    end
                    
                    initial=par.theta_z;
                    ct1=clock;
                    if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                        if sum(obj.stem_EM_options.block_tapering_block_size)==0
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                        else
                            if not(isempty(obj.stem_model.stem_data.stem_gridlist_p.tap))
                                error('Block tapering cannot be used when tapering is already enabled');
                            else
                                block_size=obj.stem_EM_options.block_tapering_block_size;
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta_block(x,par.v_z,par.correlation_type,obj.stem_model.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,block_size),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            end
                        end
                        st_par_em_step.theta_z=exp(min_result);
                    else
                        if not(isempty(obj.stem_model.stem_data.stem_gridlist_p.tap))
                            error('Block tapering cannot be used when tapering is already enabled');
                        else
                            dim=obj.stem_model.stem_data.stem_varset_p.dim;
                            dim=ones(1,par.p)*dim(1);
                            block_base=[0 cumsum(dim)];
                            theta_z_temp=st_par_em_step.theta_z;
                            if obj.stem_EM_options.workers==1
                                for b=1:length(block_base)-1
                                    idx_block=block_base(b)+1:block_base(b+1);
                                    if iscell(obj.stem_model.DistMat_z)
                                        DistMat_z_block=obj.stem_model.DistMat_z;
                                        DistMat_z_block{1}=DistMat_z_block{1}(idx_block,idx_block);
                                        DistMat_z_block{2}=DistMat_z_block{2}(idx_block,idx_block);
                                    else
                                        DistMat_z_block=obj.stem_model.DistMat_z(idx_block,idx_block);
                                    end
                                    if sum(obj.stem_EM_options.block_tapering_block_size)==0
                                        %note the input argument dim(1) because for the f-HDGM model there is a theta for each base
                                        min_result{b} = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z(b,b),par.correlation_type,DistMat_z_block,...
                                            dim(1),temp(idx_block,idx_block),T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial(:,b)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    else
                                        block_size=obj.stem_EM_options.block_tapering_block_size;
                                        %note the input argument dim(1) because for the f-HDGM model there is a theta for each base
                                        min_result{b} = fminsearch(@(x) stem_EM.geo_coreg_function_theta_block(x,par.v_z(b,b),par.correlation_type,DistMat_z_block,...
                                            dim(1),temp(idx_block,idx_block),T,block_size),log(initial(:,b)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                    theta_z_temp(:,b)=exp(min_result{b});
                                end
                            else
                                parfor b=1:length(block_base)-1
                                    idx_block=block_base(b)+1:block_base(b+1);
                                    if iscell(obj.stem_model.DistMat_z)
                                        DistMat_z_block=obj.stem_model.DistMat_z;
                                        DistMat_z_block{1}=DistMat_z_block{1}(idx_block,idx_block);
                                        DistMat_z_block{2}=DistMat_z_block{2}(idx_block,idx_block);
                                    else
                                        DistMat_z_block=obj.stem_model.DistMat_z(idx_block,idx_block);
                                    end
                                    if sum(obj.stem_EM_options.block_tapering_block_size)==0
                                        %note the input argument dim(1) because for the f-HDGM model there is a theta for each base
                                        min_result{b} = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z(b,b),par.correlation_type,DistMat_z_block,...
                                            dim(1),temp(idx_block,idx_block),T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial(:,b)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    else
                                        block_size=obj.stem_EM_options.block_tapering_block_size;
                                        %note the input argument dim(1) because for the f-HDGM model there is a theta for each base
                                        min_result{b} = fminsearch(@(x) stem_EM.geo_coreg_function_theta_block(x,par.v_z(b,b),par.correlation_type,DistMat_z_block,...
                                            dim(1),temp(idx_block,idx_block),T,block_size),log(initial(:,b)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                    theta_z_temp(:,b)=exp(min_result{b});
                                end
                            end
                            st_par_em_step.theta_z=theta_z_temp;
                        end
                    end
                    ct2=clock;
                    disp(['    theta_z updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_bp, theta_b and v_b            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_bp))
                disp('    alpha_bp updating started...');
                ct1=clock;
                alpha_bp=zeros(size(st_par_em_step.alpha_bp));
                
                for r=1:obj.stem_model.stem_par.p
                    [aj_bp_b,j_b] = obj.stem_model.get_jbp(r);
                    sum_num=0;
                    sum_den=0;
                    for t=1:T
                        if obj.stem_model.stem_data.X_bp_tv
                            tBP=t;
                        else
                            tBP=1;
                        end
                        if obj.stem_model.stem_data.X_z_tv
                            tT=t;
                        else
                            tT=1;
                        end
                        Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                        temp1=E_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp{tBP},'l'),aj_bp_b,'l');
                        temp2=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t)',M,'r'),obj.stem_model.stem_data.X_bp{tBP},'r'),j_b,'r');
                        sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                        
                        if par.p>0
                            if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                                warning('HDGM has the estimation of alpha_bp')
                                temp=obj.stem_model.stem_data.X_z{tT};
                                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                            else
                                X_z_orlated=[obj.stem_model.stem_data.X_z{tT};zeros(N-size(obj.stem_model.stem_data.X_z{tT},1),size(obj.stem_model.stem_data.X_z{tT},2))];
                            end
                            %X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                            
                            temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp{tBP},'l'),j_b,'l');
                            temp2=zeros(size(temp1,1),1);
                            if obj.stem_model.product_step>0
                                blocks=0:obj.stem_model.product_step:size(temp1,1);
                                if not(blocks(end)==size(temp1,1))
                                    blocks=cat(2,blocks,size(temp1,1));
                                end
                                for i=1:length(blocks)-1
                                    temp2(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                                end
                            else
                                temp2=diag(temp1*X_z_orlated');
                            end
                            sum_num=sum_num-sum(temp2(Lt));
                        end
                        
                        if par.k>0
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            
                            aj_p=[ones(N-Nb,1); zeros(Nb,1)];
                            
                            for k=1:K
                                temp1=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp{tBP},'l'),j_b,'l');
                                temp2=[obj.stem_model.stem_data.X_p{tP}(:,k);zeros(size(temp1,1)-size(obj.stem_model.stem_data.X_p{tP}(:,k),1),1)];
                                temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',temp2,'r'),aj_p,'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                        end
                        
                        temp1=E_wb_y1(:,t).^2+diag_Var_wb_y1(:,t);
                        temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(temp1,M,'l'),obj.stem_model.stem_data.X_bp{tBP},'b'),j_b,'b');
                        sum_den=sum_den+sum(temp1(Lt));
                    end
                    alpha_bp(r,1)=sum_num/sum_den;
                end
                
                st_par_em_step.alpha_bp=alpha_bp;
                
                ct2=clock;
                disp(['    alpha_bp updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                
                
                if obj.stem_EM_options.wpar_estimation==1
                    disp('    v_b updating started...');
                    ct1=clock;
                    
                    if Nb<=obj.stem_EM_options.mstep_system_size
                        temp=zeros(size(sum_Var_wb_y1));
                        for t=1:T
                            temp=temp+E_wb_y1(:,t)*E_wb_y1(:,t)';
                        end
                        temp=temp+sum_Var_wb_y1;
                    end
                    
                    if par.stem_par_constraints.pixel_correlated
                        %indices are permutated in order to avoid deadlock
                        v_temp=par.v_b;
                        kindex=randperm(size(par.v_b,1));
                        for k=kindex
                            hindex=randperm(size(par.v_b,1)-k)+k;
                            for h=hindex
                                initial=par.v_b(k,h);
                                if Nb<=obj.stem_EM_options.mstep_system_size
                                    min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.DistMat_b,...
                                        obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                else
                                    disp('WARNING: this operation will take a long time');
                                    min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.DistMat_b,...
                                        obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                end
                                v_temp(k,h)=min_result;
                                v_temp(h,k)=min_result;
                            end
                        end
                        if min(eig(v_temp))>=0
                            st_par_em_step.v_b=v_temp;
                        else
                            disp('    v_b is not positive definited. The last v_b is retained');
                        end
                        ct2=clock;
                        disp(['    v_b updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                    else
                        %nothing because V in this case is the identity matrix
                    end
                    
                    disp('    theta_b updating started...');
                    ct1=clock;
                    initial=par.theta_b;
                    if par.stem_par_constraints.pixel_correlated
                        if Nb<=obj.stem_EM_options.mstep_system_size
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.DistMat_b,...
                                obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            st_par_em_step.theta_b=exp(min_result);
                        else
                            if obj.stem_model.stem_data.stem_varset_b.nvar>1
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.DistMat_b,...
                                    obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                st_par_em_step.theta_b=exp(min_result);
                            else
                                s=ceil(Nb/obj.stem_EM_options.mstep_system_size);
                                step=ceil(Nb/s);
                                blocks=0:step:Nb;
                                if not(blocks(end)==Nb)
                                    blocks=[blocks Nb];
                                end
                                for j=1:length(blocks)-1
                                    block_size=blocks(j+1)-blocks(j);
                                    idx=blocks(j)+1:blocks(j+1);
                                    temp=zeros(block_size);
                                    for t=1:T
                                        temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                    end
                                    temp=temp+sum_Var_wb_y1(idx,idx);
                                    min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.DistMat_b(idx,idx),...
                                        length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                end
                                st_par_em_step.theta_b=exp(mean(min_result));
                            end
                        end
                    else
                        if Nb<=obj.stem_EM_options.mstep_system_size
                            blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                            for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),...
                                    obj.stem_model.stem_data.stem_varset_b.dim(i),temp(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                st_par_em_step.theta_b(:,i)=exp(min_result);
                            end
                        else
                            blocks_var=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                            for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                                s=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/obj.stem_EM_options.mstep_system_size);
                                step=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/s);
                                blocks=blocks_var(i):step:blocks_var(i+1);
                                if not(blocks(end)==blocks_var(i+1))
                                    blocks=cat(2,blocks,blocks_var(i+1));
                                end
                                min_result=[];
                                for j=1:length(blocks)-1
                                    block_size=blocks(j+1)-blocks(j);
                                    idx=blocks(j)+1:blocks(j+1);
                                    temp=zeros(block_size);
                                    for t=1:T
                                        temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                    end
                                    temp=temp+sum_Var_wb_y1(idx,idx);
                                    min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.DistMat_b(idx,idx),...
                                        length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                end
                                st_par_em_step.theta_b(:,i)=exp(mean(min_result));
                            end
                        end
                    end
                    ct2=clock;
                    disp(['    theta_b updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %        v_p and theta_p         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_p))
                
                if obj.stem_EM_options.wpar_estimation==1
                    disp('    v_p and theta_p updating started...');
                    ct1=clock;
                    v_temp=par.v_p;
                    
                    if size(par.v_p,3)==1&&par.k>1 %common V matrix
                        temp=zeros([size(sum_Var_wp_y1{1}),K]);
                        for z=1:K
                            for t=1:T
                                temp(:,:,z)=temp(:,:,z)+E_wp_y1(:,t,z)*E_wp_y1(:,t,z)';
                            end
                            temp(:,:,z)=temp(:,:,z)+sum_Var_wp_y1{z};
                        end
                        %indices are permutated in order to avoid deadlock
                        kindex=randperm(size(par.v_p,1));
                        for k=kindex
                            hindex=randperm(size(par.v_p,1)-k)+k;
                            for h=hindex
                                initial=par.v_p(k,h);
                                ctv1=clock;
                                if not(obj.stem_model.stem_data.stem_modeltype.is('Emulator'))
                                    min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement_shared(x,k,h,par.v_p,par.theta_p,par.correlation_type,obj.stem_model.DistMat_p,...
                                        obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                else
                                    error('Estimation of a shared V matrix is not supported when the model is of type ''Emulator''');
                                end
                                ctv2=clock;
                                disp(['    v_p(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                                v_temp(k,h)=min_result;
                                v_temp(h,k)=min_result;
                            end
                        end
                        if min(eig(v_temp))>=0
                            st_par_em_step.v_p=v_temp;
                        else
                            disp('    v_p is not positive definited. The last matrix is retained.');
                        end
                    else
                        for z=1:K
                            temp=zeros(size(sum_Var_wp_y1{z}));
                            for t=1:T
                                temp=temp+E_wp_y1(:,t,z)*E_wp_y1(:,t,z)';
                            end
                            temp=temp+sum_Var_wp_y1{z};
                            
                            %indices are permutated in order to avoid deadlock
                            kindex=randperm(size(par.v_p(:,:,z),1));
                            for k=kindex
                                hindex=randperm(size(par.v_p(:,:,z),1)-k+1)+k-1;
                                for h=hindex
                                    initial=par.v_p(k,h,z);
                                    ctv1=clock;
                                    if not(obj.stem_model.stem_data.stem_modeltype.is('Emulator'))
                                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_p(:,:,z),par.theta_p(z),par.correlation_type,obj.stem_model.DistMat_p,...
                                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    else
                                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement_emulator(x,k,h,par.v_p(:,:,z),par.theta_p,par.correlation_type,obj.stem_model.DistMat_p,...
                                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                                    end
                                    ctv2=clock;
                                    disp(['    v_p(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                                    v_temp(k,h,z)=min_result;
                                    v_temp(h,k,z)=min_result;
                                end
                            end
                        end
                        for z=1:K
                            if min(eig(v_temp(:,:,z)))>=0
                                st_par_em_step.v_p(:,:,z)=v_temp(:,:,z);
                            else
                                disp(['    v_p(:,:,',num2str(z),') is not positive definited. The last matrix is retained.']);
                            end
                        end
                    end
                    
                    if not(obj.stem_model.stem_data.stem_modeltype.is('Emulator'))
                        for z=1:K
                            temp=zeros(size(sum_Var_wp_y1{z}));
                            for t=1:T
                                temp=temp+E_wp_y1(:,t,z)*E_wp_y1(:,t,z)';
                            end
                            temp=temp+sum_Var_wp_y1{z};
                            
                            initial=par.theta_p(z);
                            ctv1=clock;
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,min(z,size(par.v_p,3))),par.correlation_type,obj.stem_model.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            st_par_em_step.theta_p(z)=exp(min_result);
                            ctv2=clock;
                            disp(['    theta_p(',num2str(z),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                        end
                    else
                        temp=zeros(size(sum_Var_wp_y1{1}));
                        for t=1:T
                            temp=temp+E_wp_y1(:,t,1)*E_wp_y1(:,t,1)';
                        end
                        temp=temp+sum_Var_wp_y1{1};
                        
                        d=length(obj.stem_model.DistMat_p);
                        for z=1:d
                            initial=par.theta_p(z);
                            ctv1=clock;
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta_emulator(x,z,par.theta_p,par.v_p(:,:,min(z,size(par.v_p,3))),par.correlation_type,obj.stem_model.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                            st_par_em_step.theta_p(z)=exp(min_result);
                            ctv2=clock;
                            disp(['    theta_p(',num2str(z),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                        end
                    end
                    ct2=clock;
                    disp(['    v_p and theta_p updating ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            
            model_changed=0;
            if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                if not(isempty(st_kalmansmoother_result))
                    clear E_e_y1;
                    clear diag_Var_e_y1;
                    clear sigma_eps;
                    clear inv_sigma_eps;
                    clear temp1;
                    clear d;
                    clear I;
                    clear K;
                    if size(obj.stem_model.stem_data.X_z{1},3)==1
                        %correlation computation
                        for i=1:N
                            L=not(isnan(obj.stem_model.stem_data.Y(i,:)));
                            a=obj.stem_model.stem_data.Y(i,L)';
                            b=st_kalmansmoother_result.zk_s(:,2:end)';
                            b=b(L,:);
                            if not(isempty(b))
                                temp=corr(a,b);
                            else
                                temp=repmat(0.0001,1,par.p);
                            end
                            if obj.stem_model.stem_data.stem_modeltype.clustering_type_is('Monopoles')
                                temp(temp<=0)=0.0001;
                            end
                            temp(isnan(temp))=0.0001;
                            obj.stem_model.stem_data.X_z{1}(i,:)=temp;
                        end
                        
                        %weight computation
                        for h=1:iteration
                            if obj.stem_model.stem_data.stem_modeltype.clustering_type_is('Monopoles')
                                obj.stem_model.stem_data.X_z{1}=obj.stem_model.stem_data.X_z{1}.^2;
                                ss=sum(obj.stem_model.stem_data.X_z{1},2);
                            else
                                obj.stem_model.stem_data.X_z{1}=obj.stem_model.stem_data.X_z{1}.^2.*sign(obj.stem_model.stem_data.X_z{1});
                                ss=sum(abs(obj.stem_model.stem_data.X_z{1}),2);
                            end
                            for j=1:size(obj.stem_model.stem_data.X_z{1},2)
                                obj.stem_model.stem_data.X_z{1}(:,j)=obj.stem_model.stem_data.X_z{1}(:,j)./ss;
                            end
                        end
                        
                        %check for empty columns
                        empty=find(sum(abs(obj.stem_model.stem_data.X_z{1}>0.01))==0);
                        if not(isempty(empty))
                            nc=st_par_em_step.p;
                            st_par_em_step.p=st_par_em_step.p-length(empty);
                            obj.stem_model.stem_data.X_z{1}(:,empty)=[];
                            sigma_eta_temp=st_par_em_step.sigma_eta;
                            sigma_eta_temp(:,empty)=[];
                            sigma_eta_temp(empty,:)=[];
                            st_par_em_step.sigma_eta=sigma_eta_temp;
                            G_temp=st_par_em_step.G;
                            G_temp(:,empty)=[];
                            G_temp(empty,:)=[];
                            st_par_em_step.G=G_temp;
                            if obj.stem_model.stem_data.stem_modeltype.is('MBC')&&obj.stem_model.stem_data.stem_modeltype.clustering_error_type_is('Dynamic')
                                sigma_eps_temp=st_par_em_step.sigma_eps;
                                sigma_eps_temp(:,empty)=[];
                                sigma_eps_temp(empty,:)=[];
                                st_par_em_step.sigma_eps=sigma_eps_temp;
                            end
                            if not(isempty(obj.stem_model.stem_data.X_beta))
                                obj.stem_model.stem_data.stem_varset_p.X_beta{1}(:,[empty,empty+nc],:)=[];
                                st_par_em_step.beta([empty,empty+nc])=[];
                            end
                            model_changed=1;
                            disp(['  Removed ',num2str(length(empty)),' empty cluster(s)']);
                        end
                    else
                        error('X_z must be a 2D matrix and not a 3D array');
                    end
                else
                    error('The Kalman smoother output is empty');
                end
                
                if obj.stem_model.stem_data.stem_modeltype.is('MBC')
                    obj.stem_model.stem_data.stem_varset_p.Y{1}=obj.stem_model.stem_data.Y;
                end
            end

            obj.stem_model.stem_par=st_par_em_step;
            ct2_mstep=clock;
            disp(['  M step ended in ',stem_misc.decode_time(etime(ct2_mstep,ct1_mstep))]);
        end
         
        function st_par_em_step = M_step_vg_and_theta(obj,E_wp_y1,sum_Var_wp_y1,index)
            %DESCRIPTION: parallel version of the M-step of the EM algorithm only for the parameters v_p and theta_p
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1) the stem_EM object
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp_k|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp_k|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp_k|Y(1)]
            %index                          - [integer >0]      (dKx1) the subset of indices from 1 to K with respect to which estimate the elements of theta_p and v_p
            %
            %OUTPUT
            %none: the stem_par property of the stem_model object is updated
            st_par_em_step=obj.stem_model.stem_par;
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            for z=index
                if Np<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wp_y1{z-index(1)+1}));
                    for t=1:size(E_wp_y1,2)
                        temp=temp+E_wp_y1(:,t,z-index(1)+1)*E_wp_y1(:,t,z-index(1)+1)';
                    end
                    temp=temp+sum_Var_wp_y1{z-index(1)+1};
                end
                kindex=randperm(size(st_par_em_step.v_p(:,:,z),1));
                for k=kindex
                    hindex=randperm(size(st_par_em_step.v_p(:,:,z),1)-k)+k;
                    for h=hindex
                        initial=st_par_em_step.v_p(k,h,z);
                        ctv1=clock;
                        if Np<=obj.stem_EM_options.mstep_system_size
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,st_par_em_step.v_p(:,:,z),st_par_em_step.theta_p(:,z),st_par_em_step.correlation_type,obj.stem_model.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                        else
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,st_par_em_step.v_p(:,:,z),st_par_em_step.theta_p(:,z),st_par_em_step.correlation_type,obj.stem_model.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                        end
                        ctv2=clock;
                        disp(['    v_p(',num2str(h),',',num2str(k),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                        st_par_em_step.v_p(k,h,z)=min_result;
                        st_par_em_step.v_p(h,k,z)=min_result;
                    end
                end
                
                initial=st_par_em_step.theta_p(:,z);
                ctv1=clock;
                if Np<=obj.stem_EM_options.mstep_system_size
                    min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.DistMat_p,...
                        obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                    st_par_em_step.theta_p(:,z)=exp(min_result);
                else
                    if obj.stem_model.stem_data.stem_varset_p.nvar>1
                        disp('WARNING: this operation will take a long time');
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                        st_par_em_step.theta_p(:,z)=exp(min_result);
                    else
                        s=ceil(Np/obj.stem_EM_options.mstep_system_size);
                        step=ceil(Np/s);
                        blocks=0:step:Np;
                        if not(blocks(end)==Np)
                            blocks=cat(2,blocks,Np);
                        end
                        for j=1:length(blocks)-1
                            block_size=blocks(j+1)-blocks(j);
                            idx=blocks(j)+1:blocks(j+1);
                            temp=zeros(block_size);
                            for t=1:size(E_wp_y1,2)
                                temp=temp+E_wp_y1(idx,t,z-index(1)+1)*E_wp_y1(idx,t,z-index(1)+1)';
                            end
                            temp=temp+sum_Var_wp_y1{z-index(1)+1}(idx,idx);
                            min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.DistMat_p(idx,idx),...
                                length(idx),temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',fminsearch_max_iter,'TolFun',0.01,'UseParallel','always'));
                        end
                        st_par_em_step.theta_p(:,z)=exp(mean(min_result));
                    end
                end
                ctv2=clock;
                disp(['    theta_p(',num2str(z),') updating ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
            end
        end
        
        %Class set function
        function set.stem_model(obj,stem_model)
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('You have to provide an object of class stem_model');
            end
        end
    end
    
    
    methods (Static)
        
        function f = geo_coreg_function_theta(log_theta,v,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to the theta_b or theta_p parameter
            %
            %INPUT
            %log_theta          - [double]      (1x1)|(2x1) natural logarithm of theta
            %v                  - [double]      (qxq) the v_b of v_q matrix
            %correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2  
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') where w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            theta=exp(log_theta);
            n_var=length(var_dims);

            if not(isempty(tapering_par))
                if not(strcmp(correlation_type,'expsphere'))
                    I=zeros(nnz(DistMat),1);
                    J=zeros(nnz(DistMat),1);
                    elements=zeros(nnz(DistMat),1);
                else
                    I=zeros(nnz(DistMat{1}),1);
                    J=zeros(nnz(DistMat{1}),1);
                    elements=zeros(nnz(DistMat{1}),1);
                end
                idx=0;
                blocks=[0 cumsum(var_dims)];
                for j=1:n_var
                    for i=j:n_var
                        if not(strcmp(correlation_type,'expsphere'))
                            B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                        else
                            for b=1:2
                                B{b} = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                            end
                        end
                        
                        corr_result=stem_misc.correlation_function(theta,B,correlation_type);
                        weights=stem_misc.wendland(B,tapering_par,correlation_type); %possibile calcolarli una sola volta???
                        corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
                        l=length(corr_result.I);
                        I(idx+1:idx+l)=corr_result.I+blocks(i);
                        J(idx+1:idx+l)=corr_result.J+blocks(j);
                        elements(idx+1:idx+l)=corr_result.correlation;
                        idx=idx+l;
                        if not(i==j)
                            I(idx+1:idx+l)=corr_result.J+blocks(j);
                            J(idx+1:idx+l)=corr_result.I+blocks(i);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                        end
                        
                    end
                end
                sigma_W=sparse(I,J,elements);
            else
                sigma_W=zeros(sum(var_dims));
                for j=1:n_var
                    for i=j:n_var
                        if not(strcmp(correlation_type,'expsphere'))
                            [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                        else
                            for b=1:2
                                [B{b},block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                            end
                        end
                        sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                        if (i~=j)
                            sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                        end
                    end
                end
            end
            if not(isempty(tapering_par))
                r = symamd(sigma_W);
                c=chol(sigma_W(r,r));
                f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
            else
                c=chol(sigma_W);
                f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
            end
        end
        
        function f = geo_coreg_function_theta_block(log_theta,v,correlation_type,DistMat,var_dims,U,T,block_size)
            %DESCRIPTION: log-likelihood evaluation with respect to the theta_b or theta_p parameter
            %
            %INPUT
            %log_theta          - [double]      (1x1) natural logarithm of theta
            %v                  - [double]      (qxq) the v_b of v_q matrix
            %correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2  
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') where w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %block_size         - [integer >0]  (1x1)|(Bx1) the size of the block-tapering block or a vector of block sizes
            %
            %OUTPUT
            %f: the log-likelihood value
            
            theta=exp(log_theta);
            n_var=length(var_dims);

            sigma_W=zeros(sum(var_dims));
            for j=1:n_var
                for i=j:n_var
                    if not(strcmp(correlation_type,'expsphere'))
                        [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                    else
                        for b=1:2
                            [B{b},block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                        end
                    end
                    sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                    if (i~=j)
                        sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                    end
                end
            end
            
            %interleaving
            if sum(sum(v-diag(diag(v))))>0
                sigma_W_int=zeros(size(sigma_W));
                U_int=zeros(size(U));
                
                n=var_dims(1);
                for i=1:n
                    for j=1:n
                        for k=1:n_var
                            for h=1:n_var
                                sigma_W_int((i-1)*n_var+k,(j-1)*n_var+h)=sigma_W((k-1)*n+i,(h-1)*n+j);
                                U_int((i-1)*n_var+k,(j-1)*n_var+h)=U((k-1)*n+i,(h-1)*n+j);
                            end
                        end
                    end
                end
            else
                sigma_W_int=sigma_W;
                U_int=U;
            end
            
            if isscalar(block_size)
                idx=0:block_size:size(sigma_W,1);
                if idx(end)<size(sigma_W,1)
                    idx=[idx size(sigma_W,1)];
                end
            else
                block_size=repmat(block_size,[1,n_var]);
                idx=[0 cumsum(block_size)];
            end
            
            f=0;
            for i=1:length(idx)-1
                c=chol(sigma_W_int(idx(i)+1:idx(i+1),idx(i)+1:idx(i+1)));
                f=f+2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U_int(idx(i)+1:idx(i+1),idx(i)+1:idx(i+1))));
            end
        end
        
        function f = geo_coreg_function_theta_emulator(log_theta,idx,theta_full,v,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to the dx1 theta_p vector when model is of the emulator type
            %
            %INPUT
            %log_theta          - [double]      (1x1) natural logarithm of the idx element of the theta vector
            %idx                - [integer>0]   (1x1) the index of the element of the theta vector which is optimized
            %theta_full         - [double]      (dx1) the full theta vector (not log transformed)
            %v                  - [double]      (qxq) the v_p matrix
            %correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      {d}x(N_p x N_p) the distance matrices
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p) sum(Var[w_p|Y(1)]+E[w_p|Y(1)]*E[w_p|Y(1)]')
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            theta=exp(log_theta);
            n_var=length(var_dims);
            
            d=length(DistMat);
            if not(isempty(tapering_par))
                I=zeros(nnz(DistMat{1}),1);
                J=zeros(nnz(DistMat{1}),1);
                elements=zeros(nnz(DistMat{1}),1);
                idx=0;
                blocks=[0 cumsum(var_dims)];
                for j=1:n_var
                    for i=j:n_var
                        B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{1});
                        weights=stem_misc.wendland(B,tapering_par,correlation_type);
                        temp=ones(length(weights),1);
                        for z=1:d
                            if z>1
                                B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{z});
                                weights=stem_misc.wendland(B,tapering_par,correlation_type);
                            end
                            if z==idx
                                corr_result=stem_misc.correlation_function(theta,B,correlation_type);
                            else
                                corr_result=stem_misc.correlation_function(theta_full(z),B,correlation_type);
                            end
                            temp=temp.*corr_result.correlation.*weights;
                            %different weights are applied for each z as the distance matrix is different for each z
                        end
                        corr_result.correlation=temp;
                        corr_result.correlation=v(i,j)*corr_result;
                        l=length(corr_result.I);
                        I(idx+1:idx+l)=corr_result.I+blocks(i);
                        J(idx+1:idx+l)=corr_result.J+blocks(j);
                        elements(idx+1:idx+l)=corr_result.correlation;
                        idx=idx+l;
                        if not(i==j)
                            I(idx+1:idx+l)=corr_result.J+blocks(j);
                            J(idx+1:idx+l)=corr_result.I+blocks(i);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                        end
                    end
                end
                sigma_W=sparse(I,J,elements);
            else
                sigma_W=ones(sum(var_dims));
                for j=1:n_var
                    for i=j:n_var
                        for z=1:d
                            [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{z});
                            if z==idx
                                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*v(i,j).*stem_misc.correlation_function(theta,B,correlation_type);
                            else
                                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*v(i,j).*stem_misc.correlation_function(theta_full(z),B,correlation_type);
                            end
                        end
                        if (i~=j)
                            sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                        end
                    end
                end
            end
            if not(isempty(tapering_par))
                r = symamd(sigma_W);
                c=chol(sigma_W(r,r));
                f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
            else
                c=chol(sigma_W);
                f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
            end
        end
        
        function f = geo_coreg_function_velement(v_element,row,col,v,theta,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %v_element          - [double]      (1x1) the extra-diagonal element of v_b or v_p
            %row                - [double]      (1x1) the row index of the v_element
            %col                - [double]      (1x1) the column index of the v_element
            %v                  - [double]      (qxq) the full v_b or v_p matrix
            %theta              - [double>0]    (1x1)|(2x1) the value of theta_b or theta_p
            %correlation type   - [string]      (1x1) correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') there w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            n_var=length(var_dims);
            v(row,col)=v_element;
            v(col,row)=v_element;
            
            if min(eig(v))>0
                if not(isempty(tapering_par))
                    sigma_W=DistMat;
                else
                    sigma_W=zeros(sum(var_dims));
                end
                
                if not(isempty(tapering_par))
                    if not(strcmp(correlation_type,'expsphere'))
                        I=zeros(nnz(DistMat),1);
                        J=zeros(nnz(DistMat),1);
                        elements=zeros(nnz(DistMat),1);
                    else
                        I=zeros(nnz(DistMat{1}),1);
                        J=zeros(nnz(DistMat{1}),1);
                        elements=zeros(nnz(DistMat{1}),1);
                    end
                    idx=0;
                    blocks=[0 cumsum(var_dims)];
                    for j=1:n_var
                        for i=j:n_var
                            if not(strcmp(correlation_type,'expsphere'))
                                B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            else
                                for b=1:2
                                    B{b} = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                                end
                            end
                            corr_result=stem_misc.correlation_function(theta,B,correlation_type);
                            weights=stem_misc.wendland(B,tapering_par,correlation_type); %possibile calcolarli una sola volta???
                            corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
                            l=length(corr_result.I);
                            I(idx+1:idx+l)=corr_result.I+blocks(i);
                            J(idx+1:idx+l)=corr_result.J+blocks(j);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                            if not(i==j)
                                I(idx+1:idx+l)=corr_result.J+blocks(j);
                                J(idx+1:idx+l)=corr_result.I+blocks(i);
                                elements(idx+1:idx+l)=corr_result.correlation;
                                idx=idx+l;
                            end
                        end
                    end
                    sigma_W=sparse(I,J,elements);
                else
                    for j=1:n_var
                        for i=j:n_var
                            if not(strcmp(correlation_type,'expsphere'))
                                [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            else
                                for b=1:2
                                    [B{b},block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                                end
                            end
                            sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                            if (i~=j)
                                sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                            end
                        end
                    end
                end
                if not(isempty(tapering_par))
                    r = symamd(sigma_W);
                    c=chol(sigma_W(r,r));
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
                else
                    c=chol(sigma_W);
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
                end
            else
                f=10^10;
            end
        end
        
        function f = geo_coreg_function_velement_block(v_element,row,col,v,theta,correlation_type,DistMat,var_dims,U,T,block_size)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %v_element          - [double]      (1x1) the extra-diagonal element of v_b or v_p
            %row                - [double]      (1x1) the row index of the v_element
            %col                - [double]      (1x1) the column index of the v_element
            %v                  - [double]      (qxq) the full v_b or v_p matrix
            %theta              - [double>0]    (1x1) the value of theta_b or theta_p
            %correlation type   - [string]      (1x1) correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') there w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %block_size         - [integer >0]  (1x1)|(Bx1) the size of the block-tapering block or a vector of block sizes
            %
            %OUTPUT
            %f: the log-likelihood value
            
            n_var=length(var_dims);
            v(row,col)=v_element;
            v(col,row)=v_element;
            
            if min(eig(v))>0
                sigma_W=zeros(sum(var_dims));
                
                for j=1:n_var
                    for i=j:n_var
                        if not(strcmp(correlation_type,'expsphere'))
                            [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                        else
                            for b=1:2
                                [B{b},block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                            end
                        end
                        sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                        if (i~=j)
                            sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                        end
                    end
                end
                
                %interleaving
                if sum(sum(v-diag(diag(v))))>0
                    sigma_W_int=zeros(size(sigma_W));
                    U_int=zeros(size(U));
                    
                    n=var_dims(1);
                    for i=1:n
                        for j=1:n
                            for k=1:n_var
                                for h=1:n_var
                                    sigma_W_int((i-1)*n_var+k,(j-1)*n_var+h)=sigma_W((k-1)*n+i,(h-1)*n+j);
                                    U_int((i-1)*n_var+k,(j-1)*n_var+h)=U((k-1)*n+i,(h-1)*n+j);
                                end
                            end
                        end
                    end
                else
                    sigma_W_int=sigma_W;
                    U_int=U;
                end
                
                if isscalar(block_size)
                    idx=0:block_size:size(sigma_W,1);
                    if idx(end)<size(sigma_W,1)
                        idx=[idx size(sigma_W,1)];
                    end
                else
                    block_size=repmat(block_size,[1,n_var]);
                    idx=[0 cumsum(block_size)];
                end
                
                f=0;
                for i=1:length(idx)-1
                    c=chol(sigma_W_int(idx(i)+1:idx(i+1),idx(i)+1:idx(i+1)));
                    f=f+2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U_int(idx(i)+1:idx(i+1),idx(i)+1:idx(i+1))));
                end
            else
                f=10^10;
            end
        end
        
        function f = geo_coreg_function_velement_shared(v_element,row,col,v,theta,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %v_element          - [double]      (1x1) the extra-diagonal element of v_p
            %row                - [double]      (1x1) the row index of the v_element
            %col                - [double]      (1x1) the column index of the v_element
            %v                  - [double]      (qxq) the full v_p matrix
            %theta              - [double>0]    (kx1) the vector theta_p
            %correlation type   - [string]      (1x1) correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      (N_p x N_p) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p x k) sum(Var[w_p|Y(1)]+E[w_p|Y(1)]*E[w_p|Y(1)]')
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            n_var=length(var_dims);
            v(row,col)=v_element;
            v(col,row)=v_element;
            K=size(U,3);
            
            if min(eig(v))>0
                f=0;
                for z=1:K
                    if not(isempty(tapering_par))
                        if not(strcmp(correlation_type,'expsphere'))
                            I=zeros(nnz(DistMat),1);
                            J=zeros(nnz(DistMat),1);
                            elements=zeros(nnz(DistMat),1);
                        else
                            I=zeros(nnz(DistMat{1}),1);
                            J=zeros(nnz(DistMat{1}),1);
                            elements=zeros(nnz(DistMat{1}),1);
                        end
                        idx=0;
                        blocks=[0 cumsum(var_dims)];
                        for j=1:n_var
                            for i=j:n_var
                                if not(strcmp(correlation_type,'expsphere'))
                                    B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                                else
                                    for b=1:2
                                        B{b} = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                                    end
                                end
                                corr_result=stem_misc.correlation_function(theta(:,z),B,correlation_type);
                                weights=stem_misc.wendland(B,tapering_par,correlation_type); 
                                corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
                                l=length(corr_result.I);
                                I(idx+1:idx+l)=corr_result.I+blocks(i);
                                J(idx+1:idx+l)=corr_result.J+blocks(j);
                                elements(idx+1:idx+l)=corr_result.correlation;
                                idx=idx+l;
                                if not(i==j)
                                    I(idx+1:idx+l)=corr_result.J+blocks(j);
                                    J(idx+1:idx+l)=corr_result.I+blocks(i);
                                    elements(idx+1:idx+l)=corr_result.correlation;
                                    idx=idx+l;
                                end
                            end
                        end
                        sigma_W=sparse(I,J,elements);
                    else
                        sigma_W=zeros(sum(var_dims));
                        for j=1:n_var
                            for i=j:n_var
                                if not(strcmp(correlation_type,'expsphere'))
                                    [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                                else
                                    for b=1:2
                                        [B{b},block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{b});
                                    end
                                end
                                sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta(z),B,correlation_type);
                                if (i~=j)
                                    sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                                end
                            end
                        end
                    end
                    if not(isempty(tapering_par))
                        r = symamd(sigma_W);
                        c=chol(sigma_W(r,r));
                        f=f+2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r,z)));
                    else
                        c=chol(sigma_W);
                        f=f+2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U(:,:,z)));
                    end
                end
            else
                f=10^10;
            end
        end
        
        function f = geo_coreg_function_velement_emulator(v_element,row,col,v,theta,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_p when the model is of the emulator type
            %
            %INPUT
            %v_element          - [double]      (1x1) the extra-diagonal element of v_p
            %row                - [double]      (1x1) the row index of the v_element
            %col                - [double]      (1x1) the column index of the v_element
            %v                  - [double]      (qxq) the full v_p matrix
            %theta              - [double>0]    (dx1) the theta_p vector
            %correlation type   - [string]      (1x1) correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      {d}x(N_p x N_p) the distance matrices
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p) sum(Var[w_p|Y(1)]+E[w_p|Y(1)]*E[w_?|Y(1)]') 
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            n_var=length(var_dims);
            v(row,col)=v_element;
            v(col,row)=v_element;
            
            if min(eig(v))>0
                if not(isempty(tapering_par))
                    sigma_W=DistMat;
                else
                    sigma_W=ones(sum(var_dims));
                end
                d=length(DistMat);
                if not(isempty(tapering_par))
                    I=zeros(nnz(DistMat{1}),1);
                    J=zeros(nnz(DistMat{1}),1);
                    elements=zeros(nnz(DistMat{1}),1);
                    idx=0;
                    blocks=[0 cumsum(var_dims)];
                    for j=1:n_var
                        for i=j:n_var
                            B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{1});
                            weights=stem_misc.wendland(B,tapering_par,correlation_type); 
                            temp=ones(length(weights),1);
                            for z=1:d
                                if z>1
                                    B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{z});
                                    weights=stem_misc.wendland(B,tapering_par,correlation_type);
                                end
                                corr_result=stem_misc.correlation_function(theta(z),B,correlation_type);
                                temp=temp.*corr_result.*weights;
                            end
                            corr_result.correlation=temp;
                            corr_result.correlation=v(i,j)*corr_result.correlation;
                            l=length(corr_result.I);
                            I(idx+1:idx+l)=corr_result.I+blocks(i);
                            J(idx+1:idx+l)=corr_result.J+blocks(j);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                            if not(i==j)
                                I(idx+1:idx+l)=corr_result.J+blocks(j);
                                J(idx+1:idx+l)=corr_result.I+blocks(i);
                                elements(idx+1:idx+l)=corr_result.correlation;
                                idx=idx+l;
                            end
                        end
                    end
                    sigma_W=sparse(I,J,elements);
                else
                    for j=1:n_var
                        for i=j:n_var
                            for z=1:d
                                [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat{z});
                                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*v(i,j).*stem_misc.correlation_function(theta(z),B,correlation_type);
                            end
                            if (i~=j)
                                sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                            end
                        end
                    end
                end
                if not(isempty(tapering_par))
                    r = symamd(sigma_W);
                    c=chol(sigma_W(r,r));
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
                else
                    c=chol(sigma_W);
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
                end
            else
                f=10^10;
            end
        end
        
        function f = function_g_element(g_element,row,col,G,sigma_eta,lambda,st_kalmansmoother_result)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %g_element                  - [double]      (1x1) the element of the diagonal of the G matrix
            %row                        - [double]      (1x1) the row index of the g_element
            %col                        - [double]      (1x1) the column index of the g_element
            %G                          - [double]      (pxp) the current full G matrix
            %sigma_eta                  - [double]      (pxp) the current full sigma_eta matrix
            %lambda                     - [double]      (1x1) the current lambda parameter
            %st_kalmansmoother_result   - [stem_kalmansmoother_result] (1x1) a stem_kalmansmoother_result object
            %
            %OUTPUT
            %f: the log-likelihood value
            
            G(row,col)=g_element;
            times=st_kalmansmoother_result.stem_datestamp.stamp;
            %set the time of the state at t_0
            m=mean(diff(times));
            times=[times(1)-m,times];
            
            zk_s=st_kalmansmoother_result.zk_s;
            Pk_s=st_kalmansmoother_result.Pk_s;
            PPk_s=st_kalmansmoother_result.PPk_s;
            
            if max(abs(eig(G)))>=1
                f=10^20;
            else
                f=0;
                for t=2:length(st_kalmansmoother_result.zk_s)
                    delta=times(t)-times(t-1);
                    n_ref=delta/lambda;
                    G_t=G.^n_ref;
                    
                    sigma_eta_t=zeros(size(sigma_eta));
                    if n_ref-round(n_ref)==0
                        for h=0:n_ref-1
                            sigma_eta_t=sigma_eta_t+G^h*sigma_eta*G^h';
                        end
                    else
                        sigma_eta_t_floor=sigma_eta_t;
                        for h=0:ceil(n_ref)-1
                            sigma_eta_t=sigma_eta_t+G^h*sigma_eta*G^h';
                            if h==floor(n_ref)-1
                                sigma_eta_t_floor=sigma_eta_t;
                            end
                        end
                        sigma_eta_t=sigma_eta_t_floor*(ceil(n_ref)-n_ref)+sigma_eta_t*(n_ref-floor(n_ref));
                    end

                    c=chol(sigma_eta_t);

                    U=(zk_s(:,t)*zk_s(:,t)'+Pk_s{t})-((zk_s(:,t)*zk_s(:,t-1)'+PPk_s{t})*G_t')-...
                      (G_t*(zk_s(:,t)*zk_s(:,t-1)'+PPk_s{t})')+(G_t*(zk_s(:,t-1)*zk_s(:,t-1)'+Pk_s{t-1})*G_t');
                    f=f+2*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
                end
            end
        end
        
        function f = function_sigma_eta_element(sigma_eta_element,row,col,G,sigma_eta,lambda,times,U)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %sigma_eta_element          - [double]      (1x1) the element of the sigma_eta matrix
            %row                        - [double]      (1x1) the row index of the g_element
            %col                        - [double]      (1x1) the column index of the g_element
            %G                          - [double]      (pxp) the current full G matrix
            %sigma_eta                  - [double]      (pxp) the current full sigma_eta matrix
            %lambda                     - [double]      (1x1) the current lambda parameter
            %times                      - [double]      (T+1x1) the times of the irregular observations
            %U                          - [double]      {T+1}(pxp) second order moments
            %
            %OUTPUT
            %f: the log-likelihood value
            
            sigma_eta(row,col)=sigma_eta_element;
            sigma_eta(col,row)=sigma_eta_element;

            av=eig(sigma_eta);
            if min(av)<=0
                f=10^20;
            else
                f=0;
                for t=2:length(times)
                    delta=times(t)-times(t-1);
                    n_ref=delta/lambda;
                    
                    sigma_eta_t=zeros(size(sigma_eta));
                    if n_ref-round(n_ref)==0
                        for h=0:n_ref-1
                            sigma_eta_t=sigma_eta_t+G^h*sigma_eta*G^h';
                        end
                    else
                        sigma_eta_t_floor=sigma_eta_t;
                        for h=0:ceil(n_ref)-1
                            sigma_eta_t=sigma_eta_t+G^h*sigma_eta*G^h';
                            if h==floor(n_ref)-1
                                sigma_eta_t_floor=sigma_eta_t;
                            end
                        end
                        sigma_eta_t=sigma_eta_t_floor*(ceil(n_ref)-n_ref)+sigma_eta_t*(n_ref-floor(n_ref));
                    end
                    
                    c=chol(sigma_eta_t);
                    f=f+2*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U{t}));
                end
            end
        end
        
    end
end

