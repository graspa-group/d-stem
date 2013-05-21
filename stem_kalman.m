%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-STEM - Distributed Space Time Expecation Maximization      %
%                                                              %
% Author: Francesco Finazzi                                    %
% E-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo - Dept. of Engineering    %
% Author website: http://www.unibg.it/pers/?francesco.finazzi  %
% Code website: https://code.google.com/p/d-stem/              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_kalman < handle
    
    %CONSTANTS
    %N   = n1_g+...+nq_g+n1_r+...+nq_r - total number of observation sites
    %N_g = n1_g+...+nq_g - total number of point sites
    %N_r = n1_r+...+nq_r - total number of pixel sites
    %N_b = n1_b+...+nq_b+n1_r+...+nq_r - total number of covariates
    %T   - number of temporal steps
    %TT = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant    
    
    properties
        stem_model=[];  %[stem_model object] (1x1) stem_model object
    end
    
    methods
        function obj = stem_kalman(stem_model)
            %DESCRIPTION: constructor of the class stem_kalman
            %
            %INPUT
            %
            %stem_model      - [stem_model object]    (1x1) stem_model object
            %
            %OUTPUT
            %obj             - [stem_kalman object]   (1x1) stem_kalman object            
            if strcmp(class(stem_model),'stem_model')
                obj.stem_model=stem_model;
            else
                error('The input argument must be of class stem_model');
            end
        end
        
        function [st_kalmanfilter_result,sigma_eps,sigma_W_r,sigma_W_g,sigma_Z,aj_rg,aj_g,M,sigma_geo] = filter(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel)
            %DESCRIPTION: Kalman filter front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]    (1x1)  stem_kalman object
            %<compute_logL>                 - [boolean]               (1x1)  (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]               (1x1)  (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]            (dTx1) (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                (1x1)  (defalut: []) full or relative path of the folder to use for distributed computation
            %    
            %OUTPUT
            %st_kalmanfilter_result         - [stem_kalmanfilter_result object] (1x1)     
            %sigma_eps                      - [double]                          (NxN) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_r                      - [double]                          (N_rxN_r) sigma_W_r matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_g                      - [double]                          {k}(N_gx_Ng) the sigma_W_g matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                          (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_rg                          - [double]                          (Nx1) the aj_rg vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_g                           - [double]                          (Nx1) the aj_g vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %M                              - [integer >0]                      (N_gx1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_geo                      - [double]                          (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            
            if nargin<2
                compute_logL=0;
            end
            if nargin<3
                enable_varcov_computation=0;
            end
            if nargin<4
                pathparallel=[];
                time_steps=[];
            end
            if nargin==4
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman filter started...');
            ct1=clock;
            
            z0=zeros(obj.stem_model.stem_par.p,1);
            P0=eye(obj.stem_model.stem_par.p);
            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;            
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;            
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter(data.Y,data.X_rg,data.X_beta,data.X_z,data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter_parallel(data.Y,data.X_rg,data.X_beta,data.X_z,data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation);
            end
            st_kalmanfilter_result = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL);
            
            ct2=clock;
            disp(['    Kalman filter ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function [st_kalmansmoother_result,sigma_eps,sigma_W_r,sigma_W_g,sigma_Z,aj_rg,aj_g,M,sigma_geo] = smoother(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel)
            %DESCRIPTION: Kalman smoother front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]    (1x1)  stem_kalman object
            %<compute_logL>                 - [boolean]               (1x1)  (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]               (1x1)  (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]            (dTx1) (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                (1x1)  (defalut: []) full or relative path of the folder to use for distributed computation
            %    
            %OUTPUT
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)     
            %sigma_eps                      - [double]                            (NxN) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_r                      - [double]                            (N_rxN_r) sigma_W_r matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_g                      - [double]                            {K}(N_gxN_g) the sigma_W_g matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                            (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_rg                          - [double]                            (Nx1) the aj_rg vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_g                           - [double]                            (Nx1) the aj_g vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %M                              - [integer >0]                        (N_gx1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_geo                      - [double]                            (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
              
            if nargin<2
              compute_logL=0;
            end
            if nargin<3
                enable_varcov_computation=0;
            end
            if nargin<4
                pathparallel=[];
                time_steps=[];
            end
            if nargin==4
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman smoother started...');
            ct1=clock;
            z0=zeros(obj.stem_model.stem_par.p,1);
            P0=eye(obj.stem_model.stem_par.p);
            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;
            [zk_s,Pk_s,PPk_s,logL] = obj.Ksmoother(data.Y,data.X_rg,data.X_beta,data.X_z,...
                data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,...
                sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,...
                time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation);
            st_kalmansmoother_result = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL,obj.stem_model.stem_data.stem_datestamp);
            ct2=clock;
            disp(['    Kalman smoother ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
    end
    
    methods (Static)
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter(Y,X_rg,X_beta,X_z,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation)
            %DESCRIPTION: Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)      the full observation matrix
            %X_rg                           - [double]     (Nx1xTT)   the full X_rg matrix
            %X_beta                         - [double]     (NxN_bxTT) the full X_beta matrix
            %X_z                         - [double]     (NxpxTT)   the full X_z matrix
            %X_g                            - [double]     (Nx1xTTxK) the full X_g matrix
            %beta                           - [double]     (N_bx1)    the beta model parameter
            %G                              - [double]     (pxp)      the G model parameter
            %sigma_eta                      - [double]     (pxp)      the sigma_eta model parameter
            %sigma_W_r                      - [double]     (N_rxN_r)  variance-covariance matrix of W_r
            %sigma_W_g                      - [double]     {K}(N_gxN_g) variance-covariance matrices of the K W_g_i
            %sigma_eps                      - [double]     (NxN)      variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)      variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_rg                          - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %aj_g                           - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_gx1)    see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)      the value of z at time t=0
            %P0                             - [double]     (pxp)      the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)      1: G and sigma_eta are diagonal matrice; 0:otherwise
            %tapering                       - [boolean]    (1x1)      1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)      1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)      1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            % 
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)    the filtered state
            %zk_u                           - [double]     (pxT+1)    the updated state
            %Pk_f                           - [double]     (pxpxT+1)  variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     (pxpxT+1)  variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)      innovation vector at time t=T
            %J                              - [double]     (pxNxT+1)  innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)      observed-data log-likelihood
            
            if nargin<20
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
           
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1));
            if enable_varcov_computation
                J_all=zeros(p,size(Y,1),T+1);
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            
            logL=0;
            for t=2:T+1
                if size(X_z,3)==1
                    tK=2;
                else
                    tK=t; %time variant
                end
                
                if compute_sigma_geo
                    if not(isempty(X_rg))
                        sigma_geo=zeros(N);
                        if size(X_rg,3)>1
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                        else
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                        end
                    end

                    if not(isempty(X_g))
                        if isempty(X_rg)
                            if tapering
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end                        
                        for k=1:size(X_g,4)
                            if size(X_g,3)>1
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                            else
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                            end
                        end
                    end
                    if isempty(X_g)&&isempty(X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                end
                
                if size(X_beta,3)==1
                    tX=2;
                else
                    tX=t; %time variant
                end
                Lt=not(isnan(Y(:,t-1))); %note the t-1
                
                X_z_orlated=X_z(:,:,tK-1);
                X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
                X_z_orlated=X_z_orlated(Lt,:);
                if stem_misc.zero_density(X_z_orlated)>90
                    X_z_orlated=sparse(X_z_orlated);
                end
                
                X_beta_orlated=X_beta(:,:,tX-1);
                X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                X_beta_orlated=X_beta_orlated(Lt,:);
                if stem_misc.zero_density(X_beta_orlated)>90
                    X_beta_orlated=sparse(X_beta_orlated);
                end                
                
                temp=sigma_geo(Lt,Lt);
                if tapering
                    if not(stem_misc.isdiagonal(temp))
                        if compute_logL
                            r = symamd(temp);
                            c=chol(temp(r,r));
                            temp2=speye(sum(Lt));
                            temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                            sigma_geo_inv=zeros(size(temp3));
                            sigma_geo_inv(r,:)=temp3;
                            clear temp2
                            clear temp3
                        end
                        temp=X_z_orlated'/temp;
                    else
                        d=1./diag(temp);
                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                        temp=X_z_orlated'*sigma_geo_inv;
                    end
                else
                    if not(stem_misc.isdiagonal(temp))
                        if compute_logL
                            c=chol(sigma_geo(Lt,Lt));
                            sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                        end
                        temp=X_z_orlated'/temp;
                    else
                        d=1./diag(temp);
                        sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                        temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                        for i=1:size(temp,1)
                            temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                        end
                    end
                end
                temp2=temp*X_z_orlated;
                
                if not(time_diagonal)
                    %filter
                    zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                    
                    %update
                    %original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    %temp=X_z_orlated'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    %temp2=temp*X_z_orlated;
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                    end
                    
                    temp3=sparse(Pk_f(:,:,t)*X_z_orlated');
                    J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                         
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1); 
                    Pk_u(:,:,t)=(eye(p)-J(:,Lt)*X_z_orlated)*Pk_f(:,:,t);  %(6.22) Stoffer
                else
                    %filter
                    zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                    
                    %update

                    %original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    %temp=X_z_orlated'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    %temp2=temp*X_z_orlated;
                    P=diag(1./diag(Pk_f(:,:,t)));
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                    end
                    temp3=Pk_f(:,:,t)*X_z_orlated';
                    J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(P+temp2)*temp;
                    
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                    Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                end
                if compute_logL
                    r = symamd(sigma_t_inv);
                    c=chol(sigma_t_inv(r,r));
                    logL=logL+(-2*sum(log(diag(c))));
                    logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                end
                clear temp
                clear temp2
                clear temp3
                if enable_varcov_computation
                    J_all(:,:,t)=J;
                end
            end
            logL=-logL/2;
            J_last=J;
            if enable_varcov_computation
                J=J_all;
            else
                J=[];
            end
        end
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter_parallel(Y,X_rg,X_beta,X_z,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)      the full observation matrix
            %X_rg                           - [double]     (Nx1xTT)   the full X_rg matrix
            %X_beta                         - [double]     (NxN_bxTT) the full X_beta matrix
            %X_z                         - [double]     (NxpxTT)   the full X_z matrix
            %X_g                            - [double]     (Nx1xTTxK) the full X_g matrix
            %beta                           - [double]     (N_bx1)    the beta model parameter
            %G                              - [double]     (pxp)      the G model parameter
            %sigma_eta                      - [double]     (pxp)      the sigma_eta model parameter
            %sigma_W_r                      - [double]     (N_rxN_r)  variance-covariance matrix of W_r
            %sigma_W_g                      - [double]     {K}(N_gxN_g) variance-covariance matrices of the K W_g_i
            %sigma_eps                      - [double]     (NxN)      variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)      variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_rg                          - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %aj_g                           - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_gx1)    see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)      the value of z at time t=0
            %P0                             - [double]     (pxp)      the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)      1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)     time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)      full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)      1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)      1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)      1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            % 
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)    the filtered state
            %zk_u                           - [double]     (pxT+1)    the updated state
            %Pk_f                           - [double]     (pxpxT+1)  variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     (pxpxT+1)  variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)      innovation vector at time t=T
            %J                              - [double]     (pxNxT+1)  innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)      observed-data log-likelihood
            
            if nargin<20
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
            
            time_steps=time_steps+1; %!!!
            
            min_ts=min(time_steps);
            max_ts=max(time_steps);
            if min_ts==2 %note the 2 due to +1
                server=1;
                %if min_ts==2 it means the call of the function is on the server
            else
                server=0;
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1));
            if enable_varcov_computation
                J_all=zeros(p,size(Y,1),T+1);
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            logL=0;
         
            if server
                for t=2:T+1
                    if size(X_z,3)==1
                        tK=2;
                    else
                        tK=t; 
                    end
                    
                    if compute_sigma_geo
                        if not(isempty(X_rg))
                            sigma_geo=zeros(N);
                            if size(X_rg,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                            end
                        end
                        
                        if not(isempty(X_g))
                            if isempty(X_rg)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_g,4)
                                if size(X_g,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_g)&&isempty(X_rg)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end
                    
                    if size(X_beta,3)==1
                        tX=2;
                    else
                        tX=t; %time variant
                    end
                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                    
                    X_z_orlated=X_z(:,:,tK-1);
                    X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
                    X_z_orlated=X_z_orlated(Lt,:);
                    if stem_misc.zero_density(X_z_orlated)>90
                        X_z_orlated=sparse(X_z_orlated);
                    end
                    
                    X_beta_orlated=X_beta(:,:,tX-1);
                    X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                    X_beta_orlated=X_beta_orlated(Lt,:);
                    if stem_misc.zero_density(X_beta_orlated)>90
                        X_beta_orlated=sparse(X_beta_orlated);
                    end
                    
                    if t>max_ts
                        %wait for the proper file from the clients
                        exit=0;
                        %disp(['        Waiting for kalman_output_',num2str(t)]);
                        while not(exit)
                            exit=exist([pathparallel,'kalman_ouput_',num2str(t),'.mat'],'file');
                        end
                        read=0;
                        while not(read)
                            try
                                load([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                read=1;
                                %disp(['        kalman_ouput_',num2str(t),' readed']);
                            catch
                            end
                            pause(0.05);
                        end
                        deleted=0;
                        while not(deleted)
                            try
                                delete([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                deleted=1;
                                %disp(['        kalman_ouput_',num2str(t),' deleted']);
                            catch
                            end
                            pause(0.05);
                        end
                    end
                    
                    if not(time_diagonal)
                        %filter
                        zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                        
                        %update
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=sigma_geo(Lt,Lt);
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                                        sigma_geo_inv=zeros(size(temp3));
                                        sigma_geo_inv(r,:)=temp3;
                                        clear temp2
                                        clear temp3
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    temp=X_z_orlated'*sigma_geo_inv;
                                end
                            else
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                                    for i=1:size(temp,1)
                                        temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                                    end
                                end
                            end
                            temp2=temp*X_z_orlated;
                            temp3=Pk_f(:,:,t)*X_z_orlated';
                            J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                        else
                            %temp and temp2 has been already reader from the file
                            temp3=Pk_f(:,:,t)*X_z_orlated';
                            J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                        end
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                        end

                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u(:,:,t)=(eye(p)-J(:,Lt)*X_z_orlated)*Pk_f(:,:,t);  %(6.22) Stoffer
                    else
                        %filter
                        zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                        
                        %update
                        
                        %original formula
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                        %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                        
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=sigma_geo(Lt,Lt);
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                                        sigma_geo_inv=zeros(size(temp3));
                                        sigma_geo_inv(r,:)=temp3;
                                        clear temp2
                                        clear temp3
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    temp=X_z_orlated'*sigma_geo_inv;
                                end
                            else
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                                    for i=1:size(temp,1)
                                        temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                                    end
                                end
                            end
                            temp2=temp*X_z_orlated;
                        else
                            %temp and temp2 has been already reader from the file
                        end
                        P=diag(1./diag(Pk_f(:,:,t)));
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                        end
                        temp3=Pk_f(:,:,t)*X_z_orlated';
                        J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(P+temp2)*temp;      
                        
                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                    end
                    if compute_logL
                        r = symamd(sigma_t_inv);
                        c=chol(sigma_t_inv(r,r));
                        logL=logL+1/(2*sum(log(diag(c))));
                        logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                    end
                    if enable_varcov_computation
                        J_all(:,:,t)=J;
                    end
                end
                logL=-logL/2;
                J_last=J;
                if enable_varcov_computation
                    J=J_all;
                else
                    J=[];
                end
            else
                %client computation
                for t=time_steps
                    if size(X_z,3)==1
                        tK=2;
                    else
                        tK=t;
                    end
                    if compute_sigma_geo
                        if not(isempty(X_rg))
                            sigma_geo=zeros(N);
                            if size(X_rg,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                            end
                        end
                        
                        if not(isempty(X_g))
                            if isempty(X_rg)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_g,4)
                                if size(X_g,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_g)&&isempty(X_rg)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end

                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                    
                    X_z_orlated=X_z(:,:,tK-1);
                    X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
                    X_z_orlated=X_z_orlated(Lt,:);
                    if stem_misc.zero_density(X_z_orlated)>90
                        X_z_orlated=sparse(X_z_orlated);
                    end
                    
                    temp=sigma_geo(Lt,Lt);
                    if tapering
                        if not(stem_misc.isdiagonal(temp))
                            temp=X_z_orlated'/temp;
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                            temp=X_z_orlated'*sigma_geo_inv;
                        end
                    else
                        if not(stem_misc.isdiagonal(temp))
                            temp=X_z_orlated'/temp;
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                            temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                            for i=1:size(temp,1)
                                temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                            end
                        end
                    end
                    temp2=temp*X_z_orlated;                    
                    
                    save([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],'temp','temp2');
                    movefile([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],[pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                    %disp(['Saved kalman_ouput_',num2str(t),'.mat']);
                end
                J_last=[];
                J=[];
            end
        end
        
        function [zk_s,Pk_s,PPk_s,logL] = Ksmoother(Y,X_rg,X_beta,X_z,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)      the full observation matrix
            %X_rg                           - [double]     (Nx1xTT)   the full X_rg matrix
            %X_beta                         - [double]     (NxN_bxTT) the full X_beta matrix
            %X_z                         - [double]     (NxpxTT)   the full X_z matrix
            %X_g                            - [double]     (Nx1xTTxK) the full X_g matrix
            %beta                           - [double]     (N_bx1)    the beta model parameter
            %G                              - [double]     (pxp)      the G model parameter
            %sigma_eta                      - [double]     (pxp)      the sigma_eta model parameter
            %sigma_W_r                      - [double]     (N_rxN_r)  variance-covariance matrix of W_r
            %sigma_W_g                      - [double]     {K}(N_gxN_g) variance-covariance matrices of the K W_g_i
            %sigma_eps                      - [double]     (NxN)      variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)      variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_rg                          - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %aj_g                           - [double]     (Nx1)      see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_gx1)    see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)      the value of z at time t=0
            %P0                             - [double]     (pxp)      the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)      1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)     time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)      full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)      1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)      1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)      1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            % 
            %OUTPUT 
            %zk_s                           - [double]     (pxT+1)    the smoothed state
            %Pk_s                           - [double]     (pxpxT+1)  variance-covariance matrix of the smoothed state
            %PPk_s                          - [double]     (pxpxT+1)  lag-one variance-covariance matrix of the smoothed state
            %logL                           - [double]     (1x1)      observed-data log-likelihood
        
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter(Y,X_rg,X_beta,X_z,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter_parallel(Y,X_rg,X_beta,X_z,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation);
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            
            H=zeros(p,p,T+1);
            Pk_s=zeros(p,p,T+1);
            zk_s=zeros(p,T+1);
            zk_s(:,end)=zk_u(:,end); %inizializzazione (6.47) Stoffer
            Pk_s(:,:,end)=Pk_u(:,:,end); %inizializzazione (6.48) Stoffer
            PPk_s=zeros(p,p,T+1);
            
            for t=T+1:-1:2
                H(:,:,t-1)=Pk_u(:,:,t-1)*G'/(Pk_f(:,:,t)); %(6.49) Stoffer
                zk_s(:,t-1)=zk_u(:,t-1)+H(:,:,t-1)*(zk_s(:,t)-zk_f(:,t)); %(6.47) Stoffer
                Pk_s(:,:,t-1)=Pk_u(:,:,t-1)+H(:,:,t-1)*(Pk_s(:,:,t)-Pk_f(:,:,t))*H(:,:,t-1)'; %(6.48) Stoffer
            end
            
            Lt=not(isnan(Y(:,end)));
            
            X_z_orlated=X_z(:,:,end);
            X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
            X_z_orlated=X_z_orlated(Lt,:);
            
            PPk_s(:,:,end)=(eye(p)-J_last(:,Lt)*X_z_orlated)*G*Pk_u(:,:,end-1); %(6.55) Stoffer
            for t=T+1:-1:3
                PPk_s(:,:,t-1)=Pk_u(:,:,t-1)*H(:,:,t-2)'+H(:,:,t-1)*(PPk_s(:,:,t)-G*Pk_u(:,:,t-1))*H(:,:,t-2)'; %(6.56) Stoffer
            end
        end
    end
end

