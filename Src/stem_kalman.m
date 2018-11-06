%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Code website: https://code.google.com/p/d-stem/                      %
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

classdef stem_kalman < handle
    
    %CONSTANTS
    %N   = n1_p+...+nq_p+n1_b+...+nq_b - total number of observation sites
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
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
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('The input argument must be of class stem_model');
            end
        end
        
        function [st_kalmanfilter_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_eta,G_tilde_diag,sigma_geo,aj_bp,aj_p,aj_z,M] = filter(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel)
            %DESCRIPTION: Kalman filter front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]              (1x1)  stem_kalman object
            %<compute_logL>                 - [boolean]                         (1x1)  (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]                         (1x1)  (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]                      (dTx1) (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                          (1x1)  (defalut: []) full or relative path of the folder to use for distributed computation
            %    
            %OUTPUT
            %st_kalmanfilter_result         - [stem_kalmanfilter_result object] (1x1)     
            %sigma_eps                      - [double]                          (NxN) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_b                      - [double]                          (N_bxN_b) sigma_W_b matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_p                      - [double]                          {k}(N_px_Ng) the sigma_W_p matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                          (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_eta                      - [double]                          (r x r) variance-covariance matrix of eta when model_name is 'HDGM' or 'f-HDGM'
            %G_tilde_diag                   - [double]                          (r x 1) diagonal of the G_tilde matrix when model_name is 'HDGM' or 'f-HDGM'
            %sigma_geo                      - [double]                          (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_bp                          - [double]                          (Nx1) the aj_bp vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_p                           - [double]                          (Nx1) the aj_p vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_z                           - [double]                          (Nx1) the aj_z vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)     
            %M                              - [integer >0]                      (N_px1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
                
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
           
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;            
            
            [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
            
            time_diagonal=obj.stem_model.stem_par.stem_par_constraints.time_diagonal;
            if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                time_diagonal=0;
            end
            tapering=obj.stem_model.tapering;  
            
            if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                s_eta=sigma_eta;
                r=length(G_tilde_diag);
                G=sparse(1:r,1:r,G_tilde_diag,r,r);
                z0=zeros(r,1);
                P0=speye(r);
            else
                s_eta=par.sigma_eta;
                G=par.G;
                z0=zeros(obj.stem_model.stem_par.p,1);
                P0=eye(obj.stem_model.stem_par.p);
            end
            
            if obj.stem_model.stem_data.stem_datestamp.irregular==1
                if not(isempty(pathparallel))
                    error('Irregular time steps not yet supported in distributed computing Kalman Filter');
                end
                times=obj.stem_model.stem_data.stem_datestamp.stamp; 
            else
                times=[];
            end
            
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter(data.Y,data.X_bp,data.X_beta,data.X_z,data.X_p,times,par.beta,G,par.lambda,s_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,obj.stem_model.stem_data.stem_modeltype);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter_parallel(data.Y,data.X_bp,data.X_beta,data.X_z,data.X_p,par.beta,G,s_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,obj.stem_model.stem_data.stem_modeltype);
            end
            st_kalmanfilter_result = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL);
            
            ct2=clock;
            disp(['    Kalman filter ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,aj_p,aj_z,M] = smoother(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel,block_tapering_block_size)
            %DESCRIPTION: Kalman smoother front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]    (1x1)         stem_kalman object
            %<compute_logL>                 - [boolean]               (1x1)         (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]               (1x1)         (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]            (dTx1)        (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                (1x1)         (defalut: []) full or relative path of the folder to use for distributed computation
            %<block_tapering_block_size>    - [integer >0]            (1x1)|(bx1)   (default: 0) the block dimension for block tapering or a vector of block sizes
            %
            %    
            %OUTPUT
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)     
            %sigma_eps                      - [double]                            (NxN|NxNxT) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_b                      - [double]                            (N_bxN_b) sigma_W_b matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_p                      - [double]                            {K}(N_pxN_p) the sigma_W_p matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                            (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_geo                      - [double]                            (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_bp                          - [double]                            (Nx1) the aj_bp vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_p                           - [double]                            (Nx1) the aj_p vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_z                           - [double]                            (Nx1) the aj_z vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details) 
            %M                              - [integer >0]                        (N_px1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            
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
            if nargin<5
                block_tapering_block_size=0;
            end
            if nargin>4
                if isempty(block_tapering_block_size)
                    block_tapering_block_size=0;
                end
            end
            if nargin==4
                error('The pathparallel input argument must be provided');
            end
            
            disp('    Kalman smoother started...');
            ct1=clock;

            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;

            time_diagonal=obj.stem_model.stem_par.stem_par_constraints.time_diagonal;
            if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                time_diagonal=0;
            end
            
            if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                s_eta=sigma_eta;
                r=length(G_tilde_diag);
                G=sparse(1:r,1:r,G_tilde_diag,r,r);
                z0=zeros(r,1);
                P0=speye(r);
            else
                s_eta=par.sigma_eta;
                G=par.G;
                z0=zeros(obj.stem_model.stem_par.p,1);
                P0=eye(obj.stem_model.stem_par.p);
            end
            
            if obj.stem_model.stem_data.stem_datestamp.irregular==1
                if not(isempty(pathparallel))
                    error('Irregular time steps not yet supported in distributed computing Kalman Smoother');
                end
                times=obj.stem_model.stem_data.stem_datestamp.stamp;
            else
                times=[];
            end
            
            if sum(block_tapering_block_size)>0
                dim=data.dim;
                if sum(diff(dim))>0
                    error('Block tapering is not yet supported in the multivariate heterotopic case');
                end
            end
            
            if sum(block_tapering_block_size)==0||not(obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                [zk_s,Pk_s,PPk_s,logL] = obj.Ksmoother(data.Y,data.X_bp,data.X_beta,data.X_z,...
                    data.X_p,times,par.beta,G,par.lambda,s_eta,sigma_W_b,...
                    sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,...
                    time_diagonal,time_steps,pathparallel,tapering,compute_logL,...
                    enable_varcov_computation,obj.stem_model.stem_data.stem_modeltype);
                st_kalmansmoother_result = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL,obj.stem_model.stem_data.stem_datestamp);
            else
                p=par.p;
                q=par.q;
                if not(p==1)
                    blocks_data=[0 cumsum(dim)];
                    if isscalar(block_tapering_block_size)
                        blocks_tapering=0:block_tapering_block_size:dim(1);
                    else
                        blocks_tapering=[0 cumsum(block_tapering_block_size)];
                    end
                    
                    if not(blocks_tapering(end)==dim(1))
                        blocks_tapering=[blocks_tapering dim(1)];
                    end
                    zk_s=zeros(r,obj.stem_model.T+1);
                    logL=0;
                    Pk_s=cell(obj.stem_model.T+1,1);
                    PPk_s=cell(obj.stem_model.T+1,1);
                    for t=1:obj.stem_model.T+1
                        Pk_s{t}=sparse(r,r);
                        PPk_s{t}=sparse(r,r);
                    end
                    sigma_Z_tap=sparse(r,r);
                    for z=1:length(blocks_tapering)-1
                        idx_tapering_q=[];
                        idx_tapering_p=[];
                        for h=1:q
                            idx_tapering_q=[idx_tapering_q (blocks_tapering(z)+1:blocks_tapering(z+1))+blocks_data(h)];
                        end
                        for h=1:p
                            idx_tapering_p=[idx_tapering_p (blocks_tapering(z)+1:blocks_tapering(z+1))+blocks_data(h)];
                        end
                        
                        Y=data.Y(idx_tapering_q,:);
                        if not(isempty(data.X_bp))
                            X_bp=data.Xbp;
                            for t=1:length(Xbp)
                                X_bp{t}=X_bp{t}(idx_tapering_q);
                            end
                        else
                            X_bp=[];
                        end
                        if not(isempty(data.X_beta))
                            X_beta=data.X_beta;
                            for t=1:length(X_beta)
                                X_beta{t}=X_beta{t}(idx_tapering_q,:);
                            end
                        else
                            X_beta=[];
                        end
                        if not(isempty(data.X_z))
                            X_z=data.X_z;
                            for t=1:length(X_z)
                                X_z{t}=X_z{t}(idx_tapering_q,idx_tapering_p);
                            end
                        else
                            X_z=[];
                        end
                        if not(isempty(data.X_p))
                            X_p=data.X_p;
                            for t=1:length(X_p)
                                X_p{t}=X_p{t}(idx_tapering_q,:);
                            end
                        else
                            X_p=[];
                        end
                        if not(isempty(G))
                            G_tap=G(idx_tapering_p,idx_tapering_p);
                        else
                            G_tap=[];
                        end
                        if not(isempty(s_eta))
                            s_eta_tap=s_eta(idx_tapering_p,idx_tapering_p);
                        else
                            s_eta_tap=[];
                        end
                        if not(isempty(sigma_W_b))
                            sigma_W_b_tap=sigma_W_b(idx_tapering_q,idx_tapering_q);
                        else
                            sigma_W_b_tap=[];
                        end
                        if not(isempty(sigma_W_p))
                            sigma_W_p_tap=sigma_W_p(idx_tapering_q,idx_tapering_q);
                        else
                            sigma_W_p_tap=[];
                        end
                        if not(isempty(sigma_eps))
                            sigma_eps_tap=sigma_eps(idx_tapering_q,idx_tapering_q);
                        else
                            sigma_eps_tap=[];
                        end
                        if not(isempty(sigma_geo))
                            sigma_geo_tap=sigma_geo(idx_tapering_q,idx_tapering_q);
                        else
                            sigma_geo_tap=[];
                        end
                        if not(isempty(aj_bp))
                            aj_bp_tap=aj_bp(idx_tapering_q);
                        else
                            aj_bp_tap=[];
                        end
                        if not(isempty(aj_p))
                            aj_p_tap=aj_p(idx_tapering_q);
                        else
                            aj_p_tap=[];
                        end
                        if not(isempty(aj_z))
                            aj_z_tap=aj_z(idx_tapering_p);
                        else
                            aj_z_tap=[];
                        end
                        if not(isempty(M))
                            M_tap=M(idx_tapering_q);
                        else
                            M_tap=[];
                        end
                        z0_tap=z0(idx_tapering_p,:);
                        P0_tap=P0(idx_tapering_p,idx_tapering_p);

                        [zk_s_tap,Pk_s_tap,PPk_s_tap,logL_tap] = obj.Ksmoother(Y,X_bp,X_beta,X_z,...
                            X_p,times,par.beta,G_tap,par.lambda,s_eta_tap,sigma_W_b_tap,...
                            sigma_W_p_tap,sigma_eps_tap,sigma_geo_tap,aj_bp_tap,aj_p_tap,aj_z_tap,M_tap,z0_tap,P0_tap,...
                            time_diagonal,time_steps,pathparallel,tapering,compute_logL,...
                            enable_varcov_computation,obj.stem_model.stem_data.stem_modeltype);

                        zk_s(idx_tapering_p,:)=zk_s_tap;
                        logL=logL+logL_tap;
                        for t=1:obj.stem_model.T+1
                            Pk_s{t}(idx_tapering_p,idx_tapering_p)=Pk_s_tap{t};
                            PPk_s{t}(idx_tapering_p,idx_tapering_p)=PPk_s_tap{t};
                        end
                        sigma_Z_tap(idx_tapering_p,idx_tapering_p)=sigma_Z(idx_tapering_p,idx_tapering_p);
                    end
                    st_kalmansmoother_result = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL,obj.stem_model.stem_data.stem_datestamp);
                    sigma_Z=sigma_Z_tap;
                else
                    error('Kalman block tapering for HDGM not yet implemented');
                end
            end
            ct2=clock;
            disp(['    Kalman smoother ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
    end
    
    methods (Static)
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter(Y,X_bp,X_beta,X_z,X_p,times,beta,G,lambda,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,stem_modeltype)
            %DESCRIPTION: Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     {T}(Nx1)    the full X_bp matrix
            %X_beta                         - [double]     {TT}(NxN_b) the full X_beta matrix
            %X_z                            - [double]     {TT}(Nxp)|{TT}(Nxr) the full X_z matrix
            %X_p                            - [double]     {TT}(NxK)   the full X_p matrix
            %times                          - [double]     (Tx1)       the vector of times of observation. Must be empty if time steps are regular
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_name is 'HDGM' or 'f-HDGM'
            %lambda                         - [double>0]   (1x1)       the lambda parameter for the case of irregular time steps
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_name is 'HDGM' or 'f-HDGM'
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN|NxNxT) variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %stem_modeltype                 - [stem_modeltype object]  (1x1)       and object of class stem_modeltype
            % 
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)                the filtered state
            %zk_u                           - [double]     (pxT+1)                the updated state
            %Pk_f                           - [double]     {T+1}(pxp)|{T+1}(rxr)  variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     {T+1}(pxp)|{T+1}(rxr)  variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)|(rxN) innovation vector at time t=T
            %J                              - [double]     {T+1}(pxN)|{T+1}(rxN)  innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)                  observed-data log-likelihood
            
            if nargin<23
                error('You have to provide all the input arguments');
            end
                        
            if not(isempty(X_beta))
                if size(X_beta{1},2)~=length(beta)
                    error('X_beta and beta are not conformable');
                end
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
            
            if not(isempty(times))
                if size(Y,2)~=length(times)
                    error(['The length of the input argument times must be ',num2str(size(Y,2))]);
                end
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
           
            if not(isempty(times))
                irregular=1;
            else
                irregular=0;
            end
            
            %set the time of the state at t_0
            if not(isempty(times))
                m=mean(diff(times));
                times=[times(1)-m,times];
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);

            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            
            Pk_f=cell(T+1,1);
            Pk_u=cell(T+1,1);
            for t=1:T+1
                if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                    Pk_f{t}=sparse(p,p);
                    Pk_u{t}=sparse(p,p);
                    J=sparse(p,size(Y,1));
                else
                    Pk_f{t}=zeros(p);
                    Pk_u{t}=zeros(p);
                    J=zeros(p,size(Y,1));
                end
            end

            if enable_varcov_computation
                J_all=cell(T+1,1);
                for t=1:T+1
                    if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                        J_all{t}=sparse(p,size(Y,1));
                    else
                        J_all{t}=zeros(p,size(Y,1));
                    end
                end
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                Pk_u{1}=sparse(P0);
            else
                Pk_u{1}=P0;
            end
            
            logL=0;
            for t=2:T+1
                if length(X_z)==1
                    tK=2;
                else
                    tK=t; %time variant
                end
                
                if compute_sigma_geo
                    if not(isempty(X_bp))
                        sigma_geo=zeros(N);
                        if length(X_bp)>1
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{t-1},'b'),aj_bp,'b');
                        else
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{1},'b'),aj_bp,'b');
                        end
                    end

                    if not(isempty(X_p))
                        if isempty(X_bp)
                            if tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end                        
                        for k=1:size(X_p{1},2)
                            if length(X_p)>1
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{t-1}(:,k),'b'),aj_p(:,k),'b');
                            else
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{1}(:,k),'b'),aj_p(:,k),'b');
                            end
                        end
                    end
                    if isempty(X_p)&&isempty(X_bp)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                end
                
                if not(isempty(X_beta))
                    if length(X_beta)==1
                        tX=2;
                    else
                        tX=t; %time variant
                    end
                end
                Lt=not(isnan(Y(:,t-1))); %note the t-1
                
                if stem_modeltype.is({'HDGM','f-HDGM'})
                    temp=X_z{tK-1};
                    X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                else
                    if N-size(X_z{1},1)>0
                        error('sparsity should be handled at this point');
                        X_z_orlated=[X_z{tK-1};zeros(N-size(X_z{tK-1},1),size(X_z{tK-1},2))];
                    else
                        X_z_orlated=X_z{tK-1};
                    end
                end
                X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');

                X_z_orlated=X_z_orlated(Lt,:);
                if stem_misc.zero_density(X_z_orlated)>90
                    X_z_orlated=sparse(X_z_orlated);
                end
                
                if not(isempty(X_beta))
                    X_beta_orlated=X_beta{tX-1};
                    X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                    X_beta_orlated=X_beta_orlated(Lt,:);
                    if stem_misc.zero_density(X_beta_orlated)>90
                        X_beta_orlated=sparse(X_beta_orlated);
                    end
                end
                
                if size(sigma_geo,3)==1
                    temp=sigma_geo(Lt,Lt);
                else
                    temp=sigma_geo(Lt,Lt,t-1);
                end
                
                if isempty(temp)
                    obs=0;
                else
                    obs=1;
                end

                if obs==1
                    if tapering
                        if not(stem_misc.isdiagonal(temp))
                            if compute_logL
                                r = symamd(temp);
                                c=chol(temp(r,r));
                                temp2=speye(sum(Lt));
                                temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
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
                                if size(sigma_geo,3)==1
                                    c=chol(sigma_geo(Lt,Lt));
                                else
                                    c=chol(sigma_geo(Lt,Lt,t-1));
                                end
                                sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                            end
                            temp=X_z_orlated'/temp;
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                            temp=X_z_orlated.*repmat(diag(sigma_geo_inv),[1,size(X_z_orlated,2)]);
                            temp=temp';
                        end
                    end
                    temp2=temp*X_z_orlated;
                end
                
                if not(time_diagonal)
                    %FILTERING
                    if irregular==1
                        delta=times(t)-times(t-1);
                        n_ref=delta/lambda;
                        G_t=G^n_ref;
                        
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

                        zk_f(:,t)=G_t*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f{t}=G_t*Pk_u{t-1}*G_t'+sigma_eta_t; %(6.20) Stoffer
                    else
                        zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f{t}=G*Pk_u{t-1}*G'+sigma_eta; %(6.20) Stoffer
                    end
                    
                    %UPDATING
                    
                    %Original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv)); 
                    %temp=X_z_orlated'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    %temp2=temp*X_z_orlated;
                    
                    if obs==1
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f{t}\eye(size(temp2)))+temp2))*temp;
                        end
                        
                        temp3=Pk_f{t}*X_z_orlated';
                        J(:,Lt)=Pk_f{t}*temp-temp3*temp'/(Pk_f{t}\speye(size(temp2))+temp2)*temp;
                        
                        if stem_misc.zero_density(J)>90
                            J=sparse(J);
                        else
                            J=full(J);
                        end
                    end
                    
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1); 
                    Pk_u{t}=(speye(p)-J(:,Lt)*X_z_orlated)*Pk_f{t};  %(6.22) Stoffer
                else
                    %FILTERING
                    if irregular==1
                        delta=times(t)-times(t-1);
                        n_ref=delta/lambda;
                        g_t=diag(G).^n_ref;
                        
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
                        
                        zk_f(:,t)=g_t.*zk_u(:,t-1); 
                        Pk_f{t}=diag((g_t.^2).*diag(Pk_u{t-1})+sigma_eta_t);
                    else
                        zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f{t}=diag(diag(G).^2.*diag(Pk_u{t-1})+diag(sigma_eta)); %(6.20) Stoffer
                    end
                    %UPDATING
                    if obs==1
                        P=spdiags(1./diag(Pk_f{t}),0,size(Pk_f{t},1),size(Pk_f{t},1));
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                        end
                        temp3=Pk_f{t}*X_z_orlated';
                        J(:,Lt)=Pk_f{t}*temp-temp3*temp'/(P+temp2)*temp;
                    end
                    
                    if stem_misc.zero_density(J)>90
                        J=sparse(J);
                    else
                        J=full(J);
                    end
                    
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                    Pk_u{t}=diag(diag((speye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f{t})); %(6.22) Stoffer
                end
                
                if compute_logL
                    if obs==1
                        if tapering
                            r = symamd(sigma_t_inv);
                            c=chol(sigma_t_inv(r,r));
                        else
                            try
                                c=chol(sigma_t_inv);
                            catch
                                a=1;
                            end
                        end
                        logL=logL+(-(2*sum(log(diag(c))))); %the negative sign is due to the fact that is the log of sigma_t that must be computed
                        logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                    end
                end
                clear temp
                clear temp2
                clear temp3
                if enable_varcov_computation
                    J_all{t}=J;
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
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter_parallel(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,stem_modeltype)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     {TT}(Nx1)   the full X_bp matrix
            %X_beta                         - [double]     {TT}(NxN_b) the full X_beta matrix
            %X_z                            - [double]     {TT}(Nxp)|{TT}(Nxr) the full X_z matrix
            %X_p                            - [double]     {TT}(NxK)   the full X_p matrix
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_name is 'HDGM' or 'f-HDGM'
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_name is 'HDGM' or 'f-HDGM'
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN|NxNxT) variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)      time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)       full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %stem_modeltype                 - [stem_modeltype object]  (1x1)    an object of class stem_modeltype
            %  
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)|(rxT+1)        the filtered state
            %zk_u                           - [double]     (pxT+1)|(rxT+1)        the updated state
            %Pk_f                           - [double]     {T+1}(pxp)|{T+1}(rxr)  variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     {T+1}(pxp)|{T+1}(rxr)  variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)|(rxN) innovation vector at time t=T
            %J                              - [double]     {T+1}(pxN)|{T+1}(pxN)  innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)                  observed-data log-likelihood
            
            if nargin<22
                error('You have to provide all the input arguments');
            end
                        
            if not(isempty(X_beta))
                if size(X_beta{1},2)~=length(beta)
                    error('X_beta and beta are not conformable');
                end
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
            
            Pk_f=cell(T+1,1);
            Pk_u=cell(T+1,1);
            for t=1:T+1
                if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                    Pk_f{t}=sparse(p,p);
                    Pk_u{t}=sparse(p,p);
                    J=sparse(p,size(Y,1));
                else
                    Pk_f{t}=zeros(p);
                    Pk_u{t}=zeros(p);
                    J=zeros(p,size(Y,1));
                end
            end

            if enable_varcov_computation
                J_all=cell(T+1,1);
                for t=1:T+1
                    if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                        J_all{t}=sparse(p,size(Y,1));
                    else
                        J_all{t}=zeros(p,size(Y,1));
                    end
                end
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                Pk_u{1}=sparse(P0);
            else
                Pk_u{1}=P0;
            end
            
            logL=0;
            if server
                for t=2:T+1
                    if length(X_z)==1
                        tK=2;
                    else
                        tK=t; 
                    end
                    
                    if compute_sigma_geo
                        if not(isempty(X_bp))
                            sigma_geo=zeros(N);
                            if length(X_bp)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{t-1},'b'),aj_bp,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{1},'b'),aj_bp,'b');
                            end
                        end
                        
                        if not(isempty(X_p))
                            if isempty(X_bp)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_p{1},4)
                                if length(X_p)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{t-1}(:,k),'b'),aj_p(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{1}(:,k),'b'),aj_p(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_p)&&isempty(X_bp)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end
                    
                    if not(isempty(X_beta))
                        if length(X_beta)==1
                            tX=2;
                        else
                            tX=t; %time variant
                        end
                    end
                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                           
                    if stem_modeltype.is({'HDGM','f-HDGM'})
                        temp=X_z{tK-1};
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[X_z{tK-1};zeros(N-size(X_z{tK-1},1),size(X_z{tK-1},2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    
                    X_z_orlated=X_z_orlated(Lt,:);
                    if stem_misc.zero_density(X_z_orlated)>90
                        X_z_orlated=sparse(X_z_orlated);
                    end
                    
                    if not(isempty(X_beta))
                        X_beta_orlated=X_beta{tX-1};
                        X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                        X_beta_orlated=X_beta_orlated(Lt,:);
                        if stem_misc.zero_density(X_beta_orlated)>90
                            X_beta_orlated=sparse(X_beta_orlated);
                        end
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
                        %FILTERING
                        zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f{t}=G*Pk_u{t-1}*G'+sigma_eta; %(6.20) Stoffer
                        
                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        %UPDATING
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            if size(sigma_geo,3)==1
                                temp=sigma_geo(Lt,Lt);
                            else
                                temp=sigma_geo(Lt,Lt,t-1);
                            end
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
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
                                        if size(sigma_geo,3)==1
                                            c=chol(sigma_geo(Lt,Lt));
                                        else
                                            c=chol(sigma_geo(Lt,Lt,t-1));
                                        end
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=X_z_orlated.*repmat(diag(sigma_geo_inv),[1,size(X_z_orlated,2)]);
                                    temp=temp';
                                end
                            end
                            temp2=temp*X_z_orlated;
                            temp3=Pk_f{t}*X_z_orlated';
                            J(:,Lt)=Pk_f{t}*temp-temp3*temp'/(Pk_f{t}\speye(size(temp2))+temp2)*temp;
                        else
                            %temp and temp2 has already been read from the file
                            temp3=Pk_f{t}*X_z_orlated';
                            J(:,Lt)=Pk_f{t}*temp-temp3*temp'/(Pk_f{t}\speye(size(temp2))+temp2)*temp;
                            if compute_logL
                                if size(sigma_geo,3)==1
                                    temp_s=sigma_geo(Lt,Lt);
                                else
                                    temp_s=sigma_geo(Lt,Lt,t-1);
                                end
                                if tapering
                                    if not(stem_misc.isdiagonal(temp))
                                        r = symamd(temp_s);
                                        c=chol(temp_s(r,r));
                                        temp2_s=speye(sum(Lt));
                                        temp3_s=full(stem_misc.chol_solve(full(c),temp2_s(r,:)));
                                        sigma_geo_inv=zeros(size(temp3_s));
                                        sigma_geo_inv(r,:)=temp3_s;
                                        clear temp2_s
                                        clear temp3_s
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                else
                                    if not(stem_misc.isdiagonal(temp))
                                        if size(sigma_geo,3)==1
                                            c=chol(sigma_geo(Lt,Lt));
                                        else
                                            c=chol(sigma_geo(Lt,Lt,t-1));
                                        end
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                end
                            end
                        end
                        
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f{t}\speye(size(temp2)))+temp2))*temp;
                            if tapering
                                r = symamd(sigma_t_inv);
                                c=chol(sigma_t_inv(r,r));
                            else
                                c=chol(sigma_t_inv);
                            end
                            logL=logL+1/(2*sum(log(diag(c))));
                            logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                        end

                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u{t}=(speye(p)-J(:,Lt)*X_z_orlated)*Pk_f{t};  %(6.22) Stoffer
                    else
                        %FILTERING
                        zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f{t}=diag(diag(G).^2.*diag(Pk_u{t-1})+diag(sigma_eta)); %(6.20) Stoffer
                        
                        %UPDATING
                        
                        %Original formula
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                        %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                        
                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            if size(sigma_geo,3)==1
                                temp=sigma_geo(Lt,Lt);
                            else
                                temp=sigma_geo(Lt,Lt,t-1);
                            end
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
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
                                        if size(sigma_geo,3)==1
                                            c=chol(sigma_geo(Lt,Lt));
                                        else
                                            c=chol(sigma_geo(Lt,Lt,t-1));
                                        end
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=X_z_orlated.*repmat(diag(sigma_geo_inv),[1,size(X_z_orlated,2)]);
                                    temp=temp';
                                end
                            end
                            temp2=temp*X_z_orlated;
                        else
                            %temp and temp2 has already been read from the file
                            if compute_logL
                                if size(sigma_geo,3)==1
                                    temp_s=sigma_geo(Lt,Lt);
                                else
                                    temp_s=sigma_geo(Lt,Lt,t-1);
                                end
                                if tapering
                                    if not(stem_misc.isdiagonal(temp))
                                        r = symamd(temp_s);
                                        c=chol(temp_s(r,r));
                                        temp2_s=speye(sum(Lt));
                                        temp3_s=full(stem_misc.chol_solve(full(c),temp2_s(r,:)));
                                        sigma_geo_inv=zeros(size(temp3_s));
                                        sigma_geo_inv(r,:)=temp3_s;
                                        clear temp2_s
                                        clear temp3_s
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                else
                                    if not(stem_misc.isdiagonal(temp))
                                        if size(sigma_geo,3)==1
                                            c=chol(sigma_geo(Lt,Lt));
                                        else
                                            c=chol(sigma_geo(Lt,Lt,t-1));
                                        end
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                end
                            end
                        end
                        
                        P=spdiags(1./diag(Pk_f{t}),0,size(Pk_f{t},1),size(Pk_f{t},1));
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                            if tapering
                                r = symamd(sigma_t_inv);
                                c=chol(sigma_t_inv(r,r));
                            else
                                c=chol(sigma_t_inv);
                            end
                            logL=logL+1/(2*sum(log(diag(c))));
                            logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                        end
                        
                        temp3=Pk_f{t}*X_z_orlated';
                        J(:,Lt)=Pk_f{t}*temp-temp3*temp'/(P+temp2)*temp;      
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u{t}=diag(diag((speye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f{t})); %(6.22) Stoffer
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
                    if length(X_z)==1
                        tK=2;
                    else
                        tK=t;
                    end
                    if compute_sigma_geo
                        if not(isempty(X_bp))
                            sigma_geo=zeros(N);
                            if length(X_bp)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{t-1},'b'),aj_bp,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp{1},'b'),aj_bp,'b');
                            end
                        end
                        
                        if not(isempty(X_p))
                            if isempty(X_bp)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_p,4)
                                if size(X_p,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{t-1}(:,k),'b'),aj_p(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p{1}(:,k),'b'),aj_p(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_p)&&isempty(X_bp)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end

                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                    
                    if stem_modeltype.is({'HDGM','f-HDGM'})
                        temp=X_z{tK-1};
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[X_z{tK-1};zeros(N-size(X_z{tK-1},1),size(X_z{tK-1},2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');

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
                            temp=X_z_orlated.*repmat(diag(sigma_geo_inv),[1,size(X_z_orlated,2)]);
                            temp=temp';
                        end
                    end
                    temp2=temp*X_z_orlated;           
                    save([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],'temp','temp2');
                    movefile([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],[pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                end
                J_last=[];
                J=[];
            end
        end
        
        function [zk_s,Pk_s,PPk_s,logL] = Ksmoother(Y,X_bp,X_beta,X_z,X_p,times,beta,G,lambda,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,stem_modeltype)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     {TT}(Nx1)   the full X_bp matrix
            %X_beta                         - [double]     {TT}(NxN_b) the full X_beta matrix
            %X_z                            - [double]     {TT}(Nxp)|{TT}(Nxr) the full X_z matrix
            %X_p                            - [double]     {TT}(NxK)   the full X_p matrix
            %times                          - [double]     (Tx1)       the vector of times of observation. Must be empty if time steps are regular
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_name is 'HDGM' or 'f-HDGM'
            %lambda                         - [double>0]   (1x1)       the lambda parameter for the case of irregular time steps
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_name is 'HDGM' or 'f-HDGM'
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN|NxNxT) variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)      time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)       full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %stem_modeltype                 - [stem_modeltype object]  (1x1)    an object of class stem_modeltype
            % 
            %OUTPUT 
            %zk_s                           - [double]     (pxT+1)|(rxT+1)       the smoothed state
            %Pk_s                           - [double]     {T+1}(pxp)|{T+1}(rxr) variance-covariance matrix of the smoothed state
            %PPk_s                          - [double]     {T+1}(pxp)|{T+1}(rxr) lag-one variance-covariance matrix of the smoothed state
            %logL                           - [double]     (1x1)                 observed-data log-likelihood
        
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter(Y,X_bp,X_beta,X_z,X_p,times,beta,G,lambda,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,stem_modeltype);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter_parallel(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,stem_modeltype);
            end
            
            if not(isempty(times))
                irregular=1;
            else
                irregular=0;
            end
            
            %set the time of the state at t_0
            if not(isempty(times))
                m=mean(diff(times));
                times=[times(1)-m,times];
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            
            H=cell(T+1,1);
            Pk_s=cell(T+1,1);
            PPk_s=cell(T+1,1);
            
            for t=1:T+1
                if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                    H{t}=sparse(p,p);
                    Pk_s{t}=sparse(p,p);
                    PPk_s{t}=sparse(p,p);
                else
                    H{t}=zeros(p);
                    Pk_s{t}=zeros(p);
                    PPk_s{t}=zeros(p);
                end
            end
            
            zk_s=zeros(p,T+1);
            zk_s(:,end)=zk_u(:,end); %inizializzazione (6.47) Stoffer
            Pk_s{end}=Pk_u{end}; %inizializzazione (6.48) Stoffer
                        
            for t=T+1:-1:2
                if irregular==1
                    delta=times(t)-times(t-1);
                    G_t=G^(delta/lambda);
                else
                    G_t=G;
                end
                H{t-1}=Pk_u{t-1}*G_t'/(Pk_f{t});
                zk_s(:,t-1)=zk_u(:,t-1)+H{t-1}*(zk_s(:,t)-zk_f(:,t)); %(6.47) Stoffer
                Pk_s{t-1}=Pk_u{t-1}+H{t-1}*(Pk_s{t}-Pk_f{t})*H{t-1}'; %(6.48) Stoffer
            end
            
            Lt=not(isnan(Y(:,end)));
            
            if stem_modeltype.is({'HDGM','f-HDGM'})&&tapering
                temp=X_z{end};
                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
            else
                if N-size(X_z,1)>0
                    X_z_orlated=[X_z{end};zeros(N-size(X_z{end},1),size(X_z{end},2))];
                else
                    X_z_orlated=X_z{end};
                end
            end
            X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
            
            X_z_orlated=X_z_orlated(Lt,:);
            if stem_misc.zero_density(X_z_orlated)>90
                X_z_orlated=sparse(X_z_orlated);
            end
            
            if irregular==1
                delta=times(t)-times(t-1);
                G_t=G^(delta/lambda);
            else
                G_t=G;
            end
            PPk_s{end}=(speye(p)-J_last(:,Lt)*X_z_orlated)*G_t*Pk_u{end-1}; %(6.55) Stoffer
            for t=T+1:-1:3
                PPk_s{t-1}=Pk_u{t-1}*H{t-1}'+H{t-1}*(PPk_s{t}-G_t*Pk_u{t-1})*H{t-2}'; %(6.56) Stoffer
            end
        end
    end
end