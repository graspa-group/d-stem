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

classdef stem_model < handle
    
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
    %
    %CONSTANTS
    %N   = n1_p+...+nq_p+n1_b+...+nq_b - total number of observation sites
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N_r - total number of elements of the z latent variable when model_name is 'HDGM'
    %r - total number of elements of the latent variable z when model_name is 'HDGM' or 'f-HDGM'
    %S = 1 if only point data are considered and S=2 if both point and pixel data are considered
    
    properties
        stem_data=[];           %[stem_data object] (1x1) object containing all the data used to estimated the model
        stem_par=[];            %[stem_par object]  (1x1) parameter set updated at each iteration of the EM algorithm
        stem_par_initial=[];    %[stem_par object]  (1x1) starting parameter set
        stem_par_sim=[];        %[stem_par object]  (1x1) parameter set used to simulate data (if data are simulated)
        estimated=0;            %[boolean] (1x1) 0: the model has not been estimated; 1: the model has been estimated
    end
    
    properties (SetAccess = private)
        DistMat_p=[];               %[double]   (N_pxN_p)|{d}x(N_pxN_p) distance matrix of the point sites
        DistMat_b=[];               %[double]   (N_bxN_b) distance matrix of the pixel sites
        DistMat_z=[];               %[double]   (N_rxN_r)|{2}x(N_pxN_p) distance matrix of the latent variable z when model_name is 'HDGM' or 'f-HDGM'. It is evaluated only if y and z are not 1:1
        
        cross_validation=0;         %[boolean]  (1x1) 0: the model has been estimated considering all the data; 1: the model has bee estimated excluding the cross-validation data.
        product_step=-1;            %[integer]  (1x1) the size of blocks when the operation diag(A*B) is performed, where A and B are matrices
        tapering=[];                %[boolean]  (1x1) 0:tapering is not enabled; 1:tapering is enabled on point sites or pixel sites
        tapering_b=[];              %[boolean]  (1x1) 0:tapering is not enabled on pixel sites; 1:tapering is enabled on pixel sites
        tapering_p=[];              %[boolean]  (1x1) 0:tapering is not enabled on point sites; 1:tapering is enabled on point sites
       
        stem_EM_result=[];          %[stem_EM_result object] (1x1) object containing all the results of the EM estimation
        stem_crossval_result=[];    %[stem_crossval_result object] {dqx1} the objects including the cross-validation results for each variable
    end
    
    methods
        
        function obj = stem_model(stem_data,stem_par)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_data - [stem_data object]  (1x1) an object of class stem_data
            %stem_par  - [stem_par object]   (1x1) an object of class stem_par
            %
            %OUTPUT
            %obj       - [stem_model object] (1x1)
            
            if nargin<2
                error('Not enough input arguments');
            end
            
            obj.stem_data=stem_data;
            obj.stem_par=stem_par;
            if not(isempty(obj.stem_data.X_beta))
                if stem_data.stem_modeltype.is('f-HDGM')&&stem_data.stem_fda.flag_beta_spline==1
                    k=getnbasis(stem_data.stem_fda.spline_basis_beta);
                    if not(size(obj.stem_data.stem_varset_p.X_beta{1},2)*k==size(obj.stem_data.X_beta{1},2))
                        error(['The number of covariates times number of basis must be equal to ',num2str(size(obj.stem_data.X_beta{1},2))]);
                    end
                else
                    if not(length(obj.stem_par.beta)==size(obj.stem_data.X_beta{1},2))
                        error(['The length of beta in stem_par must be equal to ',num2str(size(obj.stem_data.X_beta{1},2))]);
                    end
                end 
            end
            if stem_data.stem_modeltype.is('f-HDGM')
                q=length(stem_data.stem_varset_p.X_h);
                p=stem_par.p;
                X_z=cell(q,1);
                X_z_name=cell(q,1);
                for i=1:q
                    for t=1:size(stem_data.stem_varset_p.X_h{i},2)
                        temp=full(getbasismatrix(stem_data.stem_varset_p.X_h{i}(:,t),stem_data.stem_fda.spline_basis_z));
                        temp(isnan(temp))=0;
                        X_z{i}(:,:,t)=temp;
                    end
                    labels=cell(1,p);
                    for j=1:p
                        labels{1,j}=['basis_',num2str(j)];
                    end
                    X_z_name{i}=labels;
                end
                stem_data.stem_varset_p.X_z=X_z;
                stem_data.stem_varset_p.X_z_name=X_z_name;
                
                stem_data.update_data;
                 
            end
            
            if not(isempty(stem_data.stem_gridlist_b))
                if not(isempty(stem_data.stem_gridlist_b.tap))
                    obj.tapering_b=1;
                else
                    obj.tapering_b=0;
                end
            end
            if not(isempty(stem_data.stem_gridlist_p.tap))
                obj.tapering_p=1;
            else
                obj.tapering_p=0;
            end
            if not(isempty(stem_data.stem_gridlist_b))
                obj.tapering=obj.tapering_p|obj.tapering_b;
            else
                obj.tapering=obj.tapering_p;
            end
            
        end
        
        function [aj_bp] = get_aj(obj)
            %DESCRIPTION: provides the vector aj_bp used in the EM estimation
            %
            %INPUT
            %obj   - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %aj_bp - [double]            (Nx1) is the diagonal of the NxN diagonal matrix with 0 or 1 element
            %
            %NOTE
            %The elements of aj_p and aj_z from Np+1 to N are all zeros. This allows the
            %use of stem_misc.D_apply both for the pixel data and the point
            %level data avoiding the use of J_bp, J_p and J_z
            
            if not(isempty(obj.stem_data.stem_varset_b))
                
                aj_bp=zeros(obj.stem_data.N,1);
                blocks=[0 cumsum(obj.stem_data.dim)];
              
                for i=1:obj.stem_par.q
                    aj_bp(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
                for i=(obj.stem_par.q+1):obj.stem_data.nvar
                    aj_bp(blocks(i)+1:blocks(i+1))=1;
                end
            else
                aj_bp=[];
            end
            
        end
       
        function [aj_bp_b,j_b] = get_jbp(obj,r)
            %DESCRIPTION: provides the vectors aj_bp_b and j_b used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %r       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_bp_b - [double]            (Nx1) is the vector with elements equal to alpha_bp(r) only for the sites of the r-th variable
            %j_b     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            
            j_b=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
           
            for i=1:obj.stem_par.q
                if i==r
                    j_b(blocks(i)+1:blocks(i+1))=1;
                    aj_bp_b(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
            end
            for i=(obj.stem_par.q+1):obj.stem_data.nvar
                if i==r
                    j_b(blocks(i)+1:blocks(i+1))=1;
                    aj_bp_b(blocks(i)+1:blocks(i+1))=1;
                end
            end
        end
        
        function [j_p] = get_jp(obj,r)
            %DESCRIPTION: provides the vectors aj_p_bs and j_p used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %r       - [integer]           (1x1) is the index between 1 and q
            %s       - [integer]           (1x1) is the index between 1 and K
            %
            %OUTPUT
            %aj_p_bs - [double]            (Nx1) is the vector with elements equal to alpha_p(i,s) only for the sites of the r-th variable
            %j_p     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
         
            j_p=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.stem_varset_p.nvar
                if i==r
                    j_p(blocks(i)+1:blocks(i+1))=1; 
                end
            end
        end  
        
        function [j_z] = get_jz(obj,r)
            %DESCRIPTION: provides the vectors aj_z_p and j_z used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and p
            %
            %OUTPUT
            %j_z     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            if not(isempty(obj.DistMat_p))
                if not(iscell(obj.DistMat_p))
                    N=size(obj.DistMat_p,1);
                else
                    N=size(obj.DistMat_p{1},1);
                end
            end
            if not(isempty(obj.DistMat_z))
                if not(iscell(obj.DistMat_z))
                    N=size(obj.DistMat_z,1);
                else
                    N=size(obj.DistMat_z{1},1);
                end
            end
            
            j_z=zeros(N,1);
            if obj.stem_par.p==obj.stem_par.q
                blocks=[0 cumsum(obj.stem_data.dim)];
            else
                dim=obj.stem_data.dim;
                blocks=[0 cumsum(ones(1,obj.stem_par.p)*dim(1))];
            end
            for i=1:obj.stem_par.p
                if i==r
                    j_z(blocks(i)+1:blocks(i+1))=1;
                end
            end
        end  
         
        function [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,M] = get_sigma(obj)
            %DESCRIPTION: provides the variance-covariance matrices and some vectors that are used in the EM algorithm
            %
            %INPUT
            %
            %obj         - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %sigma_eps   - [double]            (NxN|NxNxT) variance-covariance matrix of epsilon
            %sigma_W_b   - [double]            (N_bxN_b) variance-covariance matrix of W_b
            %sigma_W_p   - [double]            {K}(N_pxN_p) variance-covariance matrices of the K W_p_i
            %sigma_geo   - [double]            (NxN) variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %sigma_Z     - [double]            (pxp) variance-covariance of Z
            %sigma_eta   - [double]            (r x r) variance-covariance matrix of eta when model_name is 'HDGM' or 'f-HDGM'
            %G_tilde_diag- [double]            (r x 1) diagonal of the G_tilde matrix when model_name is 'HDGM' or 'f-HDGM'
            %aj_bp       - [double]            (Nx1) see the details of the method get_jbp
            %aj_p        - [double]            (Nx1) see the details of the method get_jp
            %aj_z        - [double]            (Nx1) see the details of the method get_jz
            %M           - [integer >0]        (N_px1) see the details of the method update_M of the class stem_data
            %
            %NOTE
            %sigma_geo is provided only if it is time-invariant otherwise it is evaluated at each step of the EM algorithm

            if nargin==1
                sigma_W_b=[];
            end
            
            nvar=obj.stem_data.nvar;
            N=obj.stem_data.N;
            M=obj.stem_data.M;
            T=obj.stem_data.T;
            dim=obj.stem_data.dim;
            
            %sigma_eps
            if not(obj.stem_data.stem_modeltype.is('MBC'))
                if size(obj.stem_par.sigma_eps,3)==1
                    d=[];
                    for i=1:nvar
                        if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                            d=cat(1,d,repmat(obj.stem_par.sigma_eps(i,i),dim(i),1));
                        else
                             d=cat(1,d,repmat(obj.stem_par.sigma_eps,dim(i),1));
                        end
                    end
                    I=1:length(d);
                    sigma_eps=sparse(I,I,d);
                    
                    if obj.stem_par.flag_sigma_eps_spline==1 
                        sigma_eps = cell(T,1);
                        for t=1:T
                            if size(obj.stem_data.X_h,2)==T
                                tT=t;
                            else
                                tT=1;
                            end
                            x = full(getbasismatrix(obj.stem_data.X_h(:,tT),obj.stem_data.stem_fda.spline_basis_sigma));
                            
                            if obj.stem_par.flag_logsigma==1
                                d = exp(x*obj.stem_par.sigma_eps);
                                I=1:length(d);
                                sigma_eps{t} = sparse(I,I,d);
                            elseif obj.stem_data.stem_fda.flag_sqrsigma==1 
                                d = (x*obj.stem_par.sigma_eps).^2;
                                I=1:length(d);
                                sigma_eps{t} = sparse(I,I,d);
                            else
                                d = x*obj.stem_par.sigma_eps;
                                I=1:length(d);
                                sigma_eps{t} = sparse(I,I,d);
                            end
                           
                        end
                    end
                    
                else
                    d=[];
                    for i=1:nvar
                        d=cat(1,d,repmat(obj.stem_par.sigma_eps(i,i,:),dim(i),1));
                    end
                    sigma_eps=zeros(size(d,1),size(d,1),T);
                    for t=1:T
                        sigma_eps(:,:,t)=diag(d(:,:,t));
                    end
                end
            else
                if strcmp(obj.stem_data.stem_modeltype.clustering_error_type,'Shared')
                    d=repmat(obj.stem_par.sigma_eps(1,1),dim(1),1);
                    I=1:length(d);
                    sigma_eps=sparse(I,I,d);
                end
                if strcmp(obj.stem_data.stem_modeltype.clustering_error_type,'Proportional')
                    d=repmat(obj.stem_par.sigma_eps(1,1),dim(1),1);
                    std_variance=obj.stem_data.stem_varset_p.Y_stds{1}.^2;
                    std_variance=1./std_variance;
                    std_variance=std_variance./max(std_variance);
                    d=d.*std_variance;
                    
                    I=1:length(d);
                    sigma_eps=sparse(I,I,d);
                end
                if strcmp(obj.stem_data.stem_modeltype.clustering_error_type,'Dynamic')
                    [~,idx]=max(abs(obj.stem_data.X_z{1}),[],2);
                    d=ones(dim(1),1);
                    for i=1:obj.stem_par.p
                        L=idx==i;
                        d(L)=obj.stem_par.sigma_eps(i,i);
                    end
                    
                    I=1:length(d);
                    sigma_eps=sparse(I,I,d);
                end
            end

            %sigma_W_b
            if not(isempty(obj.stem_data.stem_varset_b))
                if not(isempty(obj.stem_data.stem_varset_b.X_bp))
                    if (nargin==1)||isempty(sigma_W_b)
                        if not(obj.tapering_b)
                            sigma_W_b=zeros(obj.stem_data.stem_varset_b.N);
                        end
                        blocks=[0 cumsum(obj.stem_data.stem_varset_b.dim)];
                        if obj.stem_par.stem_par_constraints.pixel_correlated
                            if obj.tapering_b
                                I=zeros(nnz(obj.DistMat_b),1);
                                J=zeros(nnz(obj.DistMat_b),1);
                                elements=zeros(nnz(obj.DistMat_b),1);
                            end
                            idx=0;
                            for j=1:obj.stem_data.stem_varset_b.nvar
                                for i=j:obj.stem_data.stem_varset_b.nvar
                                    idx_r=blocks(i)+1:blocks(i+1);
                                    idx_c=blocks(j)+1:blocks(j+1);
                                    if not(obj.tapering_b)
                                        sigma_W_b(idx_r,idx_c)=obj.stem_par.v_b(i,j)*stem_misc.correlation_function(...
                                            obj.stem_par.theta_b,obj.DistMat_b,obj.stem_par.correlation_type,idx_r,idx_c);
                                        if not(i==j)
                                            sigma_W_b(idx_c,idx_r)=sigma_W_b(idx_r,idx_c)';
                                        end
                                    else
                                        corr_result=stem_misc.correlation_function(obj.stem_par.theta_b,obj.DistMat_b,obj.stem_par.correlation_type,idx_r,idx_c);
                                        weights=stem_misc.wendland(obj.DistMat_b,obj.stem_data.stem_gridlist_b.tap,obj.stem_par.correlation_type,0,idx_r,idx_c);
                                        corr_result.correlation=obj.stem_par.v_b(i,j)*corr_result.correlation.*weights;
                                        siz=length(corr_result.I);
                                        I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                        J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                        elements(idx+1:idx+siz)=corr_result.correlation;
                                        idx=idx+siz;
                                        if not(i==j)
                                            I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                            J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                            elements(idx+1:idx+siz)=corr_result.correlation;
                                            idx=idx+siz;
                                        end
                                    end
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        else
                            if obj.tapering_b
                                idx=0;
                                nonzeros=0;
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    nonzeros=nonzeros+nnz(obj.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)));
                                end
                                I=zeros(nonzeros,1);
                                elements=zeros(nonzeros,1);
                            end
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                idx_rc=blocks(i)+1:blocks(i+1);
                                if not(obj.tapering_b)
                                    sigma_W_b(idx_rc,idx_rc)=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.DistMat_b,obj.stem_par.correlation_type,idx_rc,idx_rc);
                                else
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.DistMat_b,obj.stem_par.correlation_type,idx_rc,idx_rc);
                                    weights=stem_misc.wendland(obj.DistMat_b,obj.stem_data.stem_gridlist_b.tap,obj.stem_par.correlation_type,0,idx_rc,idx_rc);
                                    corr_result.correlation=obj.stem_par.v_b(i,i)*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        end
                    end
                end
            end
            clear I
            clear J
            clear elements
            clear weights
            clear corr_result
            
            %sigma_W_p
            if obj.stem_par.k>0
                blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                sigma_W_p=cell(obj.stem_par.k,1);
                if not(iscell(obj.DistMat_p))
                    for k=1:obj.stem_par.k
                        if obj.tapering_p
                            I=zeros(nnz(obj.DistMat_p),1);
                            J=zeros(nnz(obj.DistMat_p),1);
                            elements=zeros(nnz(obj.DistMat_p),1);
                        else
                            sigma_W_p{k}=obj.DistMat_p;
                        end
                        idx=0;
                        for j=1:obj.stem_data.stem_varset_p.nvar
                            for i=j:obj.stem_data.stem_varset_p.nvar
                                idx_r=blocks(i)+1:blocks(i+1);
                                idx_c=blocks(j)+1:blocks(j+1);
                                if not(obj.tapering_p)
                                    sigma_W_p{k}(idx_r,idx_c)=obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3)))*stem_misc.correlation_function(...
                                        obj.stem_par.theta_p(:,k),obj.DistMat_p,obj.stem_par.correlation_type,idx_r,idx_c);
                                    if not(i==j)
                                        sigma_W_p{k}(idx_c,idx_r)=sigma_W_p{k}(idx_r,idx_c)';
                                    end
                                else
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_p(:,k),obj.DistMat_p,obj.stem_par.correlation_type,idx_r,idx_c);
                                    weights=stem_misc.wendland(obj.DistMat_p,obj.stem_data.stem_gridlist_p.tap,obj.stem_par.correlation_type,0,idx_r,idx_c);
                                    corr_result.correlation=obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3)))*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                    if not(i==j)
                                        I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                        J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                        elements(idx+1:idx+siz)=corr_result.correlation;
                                        idx=idx+siz;
                                    end
                                end
                            end
                        end
                        if obj.tapering_p
                            sigma_W_p{k}=sparse(I,J,elements);
                        end
                    end
                else
                    d=length(obj.DistMat_p);
                    for k=1:obj.stem_par.k
                        if obj.tapering_p
                            I=zeros(nnz(obj.DistMat_p{1}),1);
                            J=zeros(nnz(obj.DistMat_p{1}),1);
                            elements=zeros(nnz(obj.DistMat_p{1}),1);
                        else
                            sigma_W_p{k}=ones(size(obj.DistMat_p{1}));
                        end
                        idx=0;
                        for j=1:obj.stem_data.stem_varset_p.nvar
                            for i=j:obj.stem_data.stem_varset_p.nvar
                                idx_r=blocks(i)+1:blocks(i+1);
                                idx_c=blocks(j)+1:blocks(j+1);
                                if not(obj.tapering_p)
                                    for z=1:d
                                        sigma_W_p{k}(idx_r,idx_c)=sigma_W_p{k}(idx_r,idx_c).*obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3))).*stem_misc.correlation_function(...
                                            obj.stem_par.theta_p(z),obj.DistMat_p{z},obj.stem_par.correlation_type,idx_r,idx_c);
                                    end
                                    if not(i==j)
                                        sigma_W_p{k}(idx_c,idx_r)=sigma_W_p{k}(idx_r,idx_c)';
                                    end
                                else
                                    weights=stem_misc.wendland(obj.DistMat_p{1},obj.stem_par.correlation_type,obj.stem_data.stem_gridlist_p.tap,0,idx_r,idx_c);
                                    temp=ones(length(weights),1);
                                    for z=1:d
                                        corr_result=stem_misc.correlation_function(obj.stem_par.theta_p(z),obj.DistMat_p{z},obj.stem_par.correlation_type,idx_r,idx_c);
                                        temp=temp.*corr_result.correlation;
                                    end
                                    corr_result.correlation=temp;
                                    
                                    corr_result.correlation=obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3)))*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                    if not(i==j)
                                        I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                        J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                        elements(idx+1:idx+siz)=corr_result.correlation;
                                        idx=idx+siz;
                                    end
                                end
                            end
                        end
                        if obj.tapering_p
                            sigma_W_p{k}=sparse(I,J,elements);
                        end
                    end
                end
            else
                sigma_W_p=[];
            end
            clear I
            clear J
            clear elements
            clear weights
            clear corr_result
            
            %sigma_eta
            if obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                if not(obj.tapering_p)
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        sigma_eta=zeros(N);
                    else
                        sigma_eta=zeros(size(obj.DistMat_z,1));
                    end
                end
                
                if obj.stem_data.stem_modeltype.is('HDGM')
                    blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                else
                    blocks=[0 cumsum(repmat(obj.stem_data.stem_varset_p.dim(1),[1,obj.stem_par.p]))];
                end
                
                if obj.tapering_p
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        I=zeros(nnz(obj.DistMat_p),1);
                        J=zeros(nnz(obj.DistMat_p),1);
                        elements=zeros(nnz(obj.DistMat_p),1);
                    else
                        I=zeros(nnz(obj.DistMat_z),1);
                        J=zeros(nnz(obj.DistMat_z),1);
                        elements=zeros(nnz(obj.DistMat_z),1);
                    end
                end
                
                idx=0;
                if not(obj.tapering_p)
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        for j=1:obj.stem_par.p
                            for i=j:obj.stem_par.p
                                idx_r=blocks(i)+1:blocks(i+1);
                                idx_c=blocks(j)+1:blocks(j+1);
                                sigma_eta(idx_r,idx_c)=obj.stem_par.v_z(i,j)*stem_misc.correlation_function(...
                                    obj.stem_par.theta_z,obj.DistMat_p,obj.stem_par.correlation_type,idx_r,idx_c);
                                if not(i==j)
                                    sigma_eta(idx_c,idx_r)=sigma_eta(idx_r,idx_c)';
                                end
                            end
                        end
                    else
                        for i=1:obj.stem_par.p
                            idx_r=blocks(i)+1:blocks(i+1);
                            idx_c=idx_r;
                            sigma_eta(idx_r,idx_c)=obj.stem_par.v_z(i,i)*stem_misc.correlation_function(...
                                obj.stem_par.theta_z(:,i),obj.DistMat_z,obj.stem_par.correlation_type,idx_r,idx_c);
                        end
                    end
                else
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        for j=1:obj.stem_par.p
                            for i=j:obj.stem_par.p
                                idx_r=blocks(i)+1:blocks(i+1);
                                idx_c=blocks(j)+1:blocks(j+1);
                                corr_result=stem_misc.correlation_function(obj.stem_par.theta_z,obj.DistMat_p,obj.stem_par.correlation_type,idx_r,idx_c);
                                weights=stem_misc.wendland(obj.DistMat_p,obj.stem_data.stem_gridlist_p.tap,obj.stem_par.correlation_type,0,idx_r,idx_c);
                            
                                corr_result.correlation=obj.stem_par.v_z(i,j)*corr_result.correlation.*weights;
                                siz=length(corr_result.I);
                                I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                elements(idx+1:idx+siz)=corr_result.correlation;
                                idx=idx+siz;
                                if not(i==j)
                                    I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                    J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                        end
                    else
                        for i=1:obj.stem_par.p
                            idx_r=blocks(i)+1:blocks(i+1);
                            idx_c=idx_r;
                            corr_result=stem_misc.correlation_function(obj.stem_par.theta_z(:,i),obj.DistMat_z,obj.stem_par.correlation_type,idx_r,idx_c);
                            weights=stem_misc.wendland(obj.DistMat_z,obj.stem_data.stem_gridlist_p.tap,obj.stem_par.correlation_type,0,idx_r,idx_c);
                            
                            corr_result.correlation=obj.stem_par.v_z(i,i)*corr_result.correlation.*weights;
                            siz=length(corr_result.I);
                            I(idx+1:idx+siz)=corr_result.I+blocks(i);
                            J(idx+1:idx+siz)=corr_result.J+blocks(i);
                            elements(idx+1:idx+siz)=corr_result.correlation;
                            idx=idx+siz;
                        end
                    end
                end
                if obj.tapering_p
                    sigma_eta=sparse(I,J,elements);
                end
            else
                sigma_eta=[];
            end
            
            %sigma_geo
            sigma_geo=[];
            aj_bp=obj.get_aj;
            if not(obj.stem_data.X_tv)
                %time invariant case
                if not(isempty(obj.stem_data.stem_varset_b))
                    sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_data.X_bp{1},'b'),aj_bp,'b');
                end
                if obj.stem_par.k>0
                    if isempty(sigma_geo)
                        %se manca il remoto allora sigma_geo non ï¿? stata
                        %ancora allocata
                        if obj.tapering
                            sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                        else
                            sigma_geo=zeros(N);
                        end
                    end
                    for k=1:obj.stem_par.k
                        sigma_geo(1:obj.stem_data.Np,1:obj.stem_data.Np)=sigma_geo(1:obj.stem_data.Np,1:obj.stem_data.Np)+...
                            stem_misc.D_apply(sigma_W_p{k},[obj.stem_data.X_p{1}(:,k);...
                            zeros(obj.stem_data.Np-size(obj.stem_data.X_p{1}(:,k),1),1)],'b');
                    end
                end
                if isempty(sigma_geo)
                    sigma_geo=sigma_eps;
                else
                    sigma_geo=sigma_geo+sigma_eps;
                end
            end
            
            if obj.stem_par.p>0
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    temp=abs(obj.stem_par.G-eye(size(obj.stem_par.G)));
                    if sum(temp(:))==0
                        %random walk
                        sigma_Z=T*obj.stem_par.sigma_eta;
                    else
                        sigma_Z=(eye(obj.stem_par.p^2)-kron(obj.stem_par.G,obj.stem_par.G))\obj.stem_par.sigma_eta(:);%variance of a VAR (Lutkepohl pag.27)
                        sigma_Z=reshape(sigma_Z,obj.stem_par.p,obj.stem_par.p);
                    end
                    G_tilde_diag=[];
                else
                    temp=abs(obj.stem_par.G-eye(size(obj.stem_par.G)));
                    if sum(temp(:))==0
                        %random walk
                        sigma_Z=T*obj.stem_par.sigma_eta;
                    else
                        dim=obj.stem_data.dim;
                        if obj.stem_data.stem_modeltype.is('f-HDGM')
                            G_tilde_diag=kron(diag(obj.stem_par.G),ones(dim(1),1));
                        else
                            G_tilde_diag=[];
                            for i=1:obj.stem_par.p
                                G_tilde_diag=cat(1,G_tilde_diag,obj.stem_par.G(i,i)*ones(dim(i),1));
                            end
                        end
                        sigma_Z=reshape(1./(1-kron(G_tilde_diag,G_tilde_diag)).*sigma_eta(:),length(G_tilde_diag),length(G_tilde_diag));
                    end
                end
            else
                sigma_Z=[];
                G_tilde_diag=[];
            end
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            %DESCRIPTION: is the front-end method to simulate data from this stem_model object
            %
            %INPUT
            %
            %obj               - [stem_model object] (1x1) the stem_model object
            %nan_rate          - [double [0,1]]      (Sx1) (default: []) missing rate of the simulated data
            %nan_pattern_par   - [double >0]         (Sx1) (default: []) (UoM: km) parameter of the exponential spatial correlation function used to define the spatial pattern of the missing data
            %
            %OUTPUT
            %
            %none: an object of class stem_sim is created. The method stem_sim.simulate actually simulates the data and it fills the stem_data object
            st_sim=stem_sim(obj);
            if nargin<2
                nan_rate=[];
                nan_pattern_par=[];
            end
            if nargin==2
                error('nan_pattern_par must also be provided');
            end
            st_sim.simulate(nan_rate,nan_pattern_par);
            obj.stem_data.simulated=1;
        end        
        
        function EM_estimate(obj,stem_EM_options)
            %DESCRIPTION: front-end method for EM estimation of this model
            %
            %INPUT
            %
            %obj                - [stem_model object]      (1x1) the stem_model object
            %stem_EM_options    - [stem_EM_options object] (1x1) an object of class stem_EM_options
            %
            %OUTPUT
            %
            %none: the stem_EM object is created. The method stem_EM.estimate is used to estimate the model and it updates the stem_par object

            obj.set_distance_matrix;
            
            if length(stem_EM_options.block_tapering_block_size)>1
                dim=obj.dim;
                if not(sum(stem_EM_options.block_tapering_block_size)==dim(1))
                    error('The sum of the elements of block_tapering_block_size must be equal to the number of spatial sites');
                end
            end

            if sum(stem_EM_options.block_tapering_block_size)>0&&not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                error('block_tapering_block_size must be 0 when modeltype is not f-HDGM')
            end
            
            if not(isempty(obj.stem_data.stem_crossval))
                if length(stem_EM_options.block_tapering_block_size)>1
                    idx=[];
                    for z=1:length(stem_EM_options.block_tapering_block_size)
                        idx=[idx ones(1,stem_EM_options.block_tapering_block_size(z))*z];
                    end
                    if obj.stem_par.stem_modeltype.is('f-HDGM')
                        idx(obj.stem_data.stem_crossval.indices)=[];
                    else
                        idx(obj.stem_data.stem_crossval.indices{1})=[];
                    end
                    stem_EM_options.block_tapering_block_size=diff([0 find(diff([idx 0]))]);
                end
                disp('Data modification for cross-validation started...');
                
                disp('******************')
                disp('Cropping sites...');
                
                if obj.stem_par.stem_modeltype.is('f-HDGM')
                   
                    obj.stem_crossval_result=stem_crossval_result();

                    idx_var=obj.stem_data.stem_varset_p.get_Y_index(obj.stem_data.stem_crossval.variable_name);
                    if isempty(idx_var)
                        error(['Cross-validation variable',obj.stem_data.stem_crossval.variable_name,'not found.']);
                    end

                    %recover the indices of the cross-validation sites
                    indices=obj.stem_data.stem_crossval.indices;
                    for idx = idx_var'
                        Y={obj.stem_data.stem_varset_p.Y{idx}(indices,:)};
                        Y_name=obj.stem_data.stem_varset_p.Y_name(idx);
                        
                        X_bp={};
                        X_bp_name={};
                        
                        if not(isempty(obj.stem_data.stem_varset_p.X_beta))
                            X_beta={obj.stem_data.stem_varset_p.X_beta{idx}(indices,:,:)};
                            X_beta_name=obj.stem_data.stem_varset_p.X_beta_name(idx);
                        else
                            X_beta={};
                            X_beta_name={};
                        end
                        if not(isempty(obj.stem_data.stem_varset_p.X_h))
                            X_h={obj.stem_data.stem_varset_p.X_h{idx}(indices,:,:)};
                            X_h_name=obj.stem_data.stem_varset_p.X_h_name;
                        else
                            X_h={};
                            X_h_name={};
                        end
                        if not(isempty(obj.stem_data.stem_varset_p.X_z))
                            X_z={obj.stem_data.stem_varset_p.X_z{idx}(indices,:,:)};
                            X_z_name=obj.stem_data.stem_varset_p.X_z_name(idx);
                        else
                            X_z={};
                            X_z_name={};
                        end
                        
                        X_p={};
                        X_p_name={};

                        obj.stem_data.stem_crossval.stem_varset{idx}=stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name,X_h,X_h_name);
                        obj.stem_data.stem_crossval.stem_gridlist{idx}=stem_gridlist();
                        coordinate=obj.stem_data.stem_gridlist_p.grid{idx}.coordinate;
                        st_grid=stem_grid(coordinate(indices,:),'deg','sparse','point');
                        obj.stem_data.stem_crossval.stem_gridlist{idx}.add(st_grid); 
                    end
                    obj.stem_data.site_crop(obj.stem_data.stem_crossval.type,obj.stem_data.stem_crossval.variable_name,indices,1,0);
                    
                else
                    for i=1:length(obj.stem_data.stem_crossval.variable_name)
                    
                        obj.stem_crossval_result{i}=stem_crossval_result();

                        idx_var=obj.stem_data.stem_varset_p.get_Y_index(obj.stem_data.stem_crossval.variable_name{i});
                        if isempty(idx_var)
                            error(['Cross-validation variable',obj.stem_data.stem_crossval.variable_name{i},'not found.']);
                        end

                        %recover the indices of the cross-validation sites
                        indices=obj.stem_data.stem_crossval.indices{i};
                        if isempty(indices)
                            Y={};
                            Y_name={};
                            X_bp={};
                            X_bp_name={};
                            
                            X_beta={};
                            X_beta_name={};
                            
                            X_h={};
                            X_h_name={};
                            X_z={};
                            X_z_name={};
                            X_p={};
                            X_p_name={};

                        else
                            Y={obj.stem_data.stem_varset_p.Y{idx_var}(indices,:)};
                            Y_name=obj.stem_data.stem_varset_p.Y_name(idx_var);
                            if not(isempty(obj.stem_data.stem_varset_p.X_bp))
                                X_bp={obj.stem_data.stem_varset_p.X_bp{idx_var}(indices,:,:)};
                                X_bp_name=obj.stem_data.stem_varset_p.X_bp(idx_var);
                            else
                                X_bp={};
                                X_bp_name={};
                            end
                            if not(isempty(obj.stem_data.stem_varset_p.X_beta))
                                X_beta={obj.stem_data.stem_varset_p.X_beta{idx_var}(indices,:,:)};
                                X_beta_name=obj.stem_data.stem_varset_p.X_beta_name(idx_var);
                            else
                                X_beta={};
                                X_beta_name={};
                            end
                            if not(isempty(obj.stem_data.stem_varset_p.X_h))
                                X_h={obj.stem_data.stem_varset_p.X_h{idx_var}(indices,:,:)};
                                X_h_name=obj.stem_data.stem_varset_p.X_h_name;
                            else
                                X_h={};
                                X_h_name={};
                            end
                            if not(isempty(obj.stem_data.stem_varset_p.X_z))
                                X_z={obj.stem_data.stem_varset_p.X_z{idx_var}(indices,:,:)};
                                X_z_name=obj.stem_data.stem_varset_p.X_z_name(idx_var);
                            else
                                X_z={};
                                X_z_name={};
                            end
                            if not(isempty(obj.stem_data.stem_varset_p.X_p))
                                X_p={obj.stem_data.stem_varset_p.X_p{idx_var}(indices,:,:,:)};
                                X_p_name=obj.stem_data.stem_varset_p.X_p_name(idx_var);
                            else
                                X_p={};
                                X_p_name={};
                            end

                            %set the cross_mindistance vector
                            if not(isempty(obj.DistMat_p))
                                dim=obj.stem_data.dim;
                                blocks=[0 cumsum(dim)];
                                if not(strcmp(obj.stem_par.correlation_type,'expsphere'))
                                    temp_dist=obj.DistMat_p(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                                else
                                    temp_dist=obj.DistMat_p{1}(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                                end
                                temp_dist=temp_dist(indices,:);
                                temp_dist(:,indices)=[];
                                obj.stem_crossval_result{i}.min_distance=min(temp_dist,[],2);
                                clear temp_dist
                            end
                        end
                        
                        obj.stem_data.stem_crossval.stem_varset{i}=stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name,X_h,X_h_name);
                        obj.stem_data.stem_crossval.stem_gridlist{i}=stem_gridlist();
                        coordinate=obj.stem_data.stem_gridlist_p.grid{idx_var}.coordinate;
                        st_grid=stem_grid(coordinate(indices,:),'deg','sparse','point');
                        obj.stem_data.stem_crossval.stem_gridlist{i}.add(st_grid);
                        %remove the cross-validation data from the estimation dataset
                        obj.stem_data.site_crop(obj.stem_data.stem_crossval.type{i},obj.stem_data.stem_crossval.variable_name{i},indices,1,0);
                    end
                end

                disp('Cropping ended.');
                disp('******************')
                obj.stem_data.update_data;
                obj.set_distance_matrix;
                obj.cross_validation=1;
                disp('Data modification ended.');
            else
                obj.cross_validation=0;
            end
            
            st_EM=stem_EM(obj,stem_EM_options);
            %set the current parameter value with the estimated initial value
            obj.stem_par=obj.stem_par_initial;
            obj.data_summary;
            
            if stem_EM_options.workers==1
                obj.stem_EM_result=st_EM.estimate();
            else
                delete(gcp('nocreate'))
                poolobj = parpool(stem_EM_options.workers);
                obj.stem_EM_result=st_EM.estimate();
                delete(poolobj);
            end
            
            obj.estimated=1;
            if obj.cross_validation
                obj.fill_crosval_result();
            end
        end
          
        function fill_crosval_result(obj)
            %DESCRIPTION: obtain cross-validation result
            %
            %INPUT
            %obj  - [stem_model object]   (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: stem_crossval_result is created in the stem_crossval_result property of the stem_model object
            
            if obj.cross_validation
                %check if the cross-validation sites are the same for all the variables
                equal=1;
                for i=1:length(obj.stem_data.stem_crossval.stem_gridlist)-1
                    if not(size(obj.stem_data.stem_crossval.stem_gridlist{i}.grid{1}.coordinate,1)==size(obj.stem_data.stem_crossval.stem_gridlist{i+1}.grid{1}.coordinate,1))
                        equal=0;
                    else
                        if sum(sum(obj.stem_data.stem_crossval.stem_gridlist{i}.grid{1}.coordinate-obj.stem_data.stem_crossval.stem_gridlist{i+1}.grid{1}.coordinate))>0
                            equal=0;
                        end
                    end
                end
                
                st_krig_options=stem_krig_options();
                st_krig_options.block_size=300;
                st_krig_options.nn_size=obj.stem_data.stem_crossval.nn_size;
                st_krig_options.back_transform=0;
                st_krig_options.no_varcov=0;
                st_krig_options.crossval=1;
                
                if equal
                    %if they are the same, only one multivariate-kriging is enough
                    disp('Kriging over cross-validation sites');
                    
                    st_krig_data=stem_krig_data(obj.stem_data.stem_crossval.stem_gridlist{1}.grid{1});
                    
                    st_krig=stem_krig(obj,st_krig_data);
                    st_krig_result=st_krig.kriging(st_krig_options);
                else
                    %if they are different, kriging is repeated with respect to the different coordinates
                    st_krig_result=cell(length(obj.stem_data.stem_crossval.variable_name),1);
                    for i=1:length(obj.stem_data.stem_crossval.variable_name)
                        if obj.stem_data.stem_crossval.stem_varset{i}.N>0
                            disp(['Kriging over cross-validation sites of variable ',obj.stem_data.stem_crossval.variable_name{i}]);
                        
                            st_krig_data=stem_krig_data(obj.stem_data.stem_crossval.stem_gridlist{i}.grid{1});

                            st_krig=stem_krig(obj,st_krig_data);
                            st_krig_result_temp=st_krig.kriging(st_krig_options,i);
                            st_krig_result{i}=st_krig_result_temp{i};
                        end
                    end
                end  
                
                if obj.stem_par.stem_modeltype.is('f-HDGM')
                    
                    obj.stem_crossval_result.res=cell(length(st_krig_result),1);
                    obj.stem_crossval_result.res_back=cell(length(st_krig_result),1);
                    obj.stem_crossval_result.y_back=cell(length(st_krig_result),1);
                    obj.stem_crossval_result.y_hat_back=cell(length(st_krig_result),1);

                    for i=1:length(st_krig_result)
                        obj.stem_crossval_result.res{i}=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}-st_krig_result{i}.y_hat;

                        if obj.stem_data.stem_varset_p.standardized
                            s=obj.stem_data.stem_varset_p.Y_stds{i};
                            m=obj.stem_data.stem_varset_p.Y_means{i};
                        else
                            s=1;
                            m=0;
                        end

                        if not(obj.stem_data.stem_varset_p.log_transformed)&&not(obj.stem_data.stem_varset_p.boxcox_transformed)
                            y_hat_back=st_krig_result{i}.y_hat*s+m;
                            y=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m;
                        end

                        if obj.stem_data.stem_varset_p.log_transformed
                            y_hat_back=st_krig_result{i}.y_hat;
                            var_y_hat=st_krig_result{i}.diag_Var_y_hat;
                            y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                            y=exp(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m);
                        end

                        if obj.stem_data.stem_varset_p.boxcox_transformed
                            lambda=obj.stem_data.stem_varset_p.Y_lambda{idx_var};

                            y_hat_back=st_krig_result{i}.y_hat;
                            var_y_hat=st_krig_result{i}.diag_Var_y_hat;

                            y_hat_back=(max(0,(lambda*(y_hat_back*s+m)+1)).^(1/lambda)).*(1+(var_y_hat*s^2*(1-lambda))./(2*(lambda*(y_hat_back*s+m)+1).^2));
                            y=(lambda*(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m)+1).^(1/lambda);
                        end

                        obj.stem_crossval_result.res_back{i}=y-y_hat_back;
                        obj.stem_crossval_result.y_back{i}=y;
                        obj.stem_crossval_result.y_hat_back{i}=y_hat_back;
                    end
                   
                    %w.r.t t
                    res=nan(obj.stem_par.q,size(obj.stem_crossval_result.res{1},1),obj.stem_data.T);
                    for qq=1:obj.stem_par.q
                        res(qq,:,:)=obj.stem_crossval_result.res{qq};
                    end
                    
                    Y=nan(obj.stem_par.q,size(obj.stem_crossval_result.res{1},1),obj.stem_data.T);
                    for qq=1:obj.stem_par.q
                        Y(qq,:,:)=obj.stem_crossval_result.y_back{qq};
                    end
                    
                    if size(res,3)==1
                        res_t=res(:);
                        Y_t=Y(:);
                    else
                        res_t=squeeze(reshape(res,1,size(res,1)*size(res,2),size(res,3)));
                        Y_t=squeeze(reshape(Y,1,size(res,1)*size(res,2),size(res,3)));
                    end
                    obj.stem_crossval_result.cv_mse_t=nanvar(res_t,1,1)';
                    relative_cv_mse_t=obj.stem_crossval_result.cv_mse_t./nanvar(Y_t,1,1)';
                    obj.stem_crossval_result.cv_R2_t=1-relative_cv_mse_t;
                    
                    %w.r.t s
                    if size(res,2)==1
                        res_s=res(:);
                        Y_s=Y(:);
                    else
                        res1=shiftdim(res,2);
                        Y1=shiftdim(Y,2);
                        res_s=squeeze(reshape(res1,1,size(res1,1)*size(res1,2),size(res1,3)));
                        Y_s=squeeze(reshape(Y1,1,size(res1,1)*size(res1,2),size(res1,3)));
                    end
                    obj.stem_crossval_result.cv_mse_s=nanvar(res_s,1,1)';
                    relative_cv_mse_s=obj.stem_crossval_result.cv_mse_s./nanvar(Y_s,1,1)';
                    obj.stem_crossval_result.cv_R2_s=1-relative_cv_mse_s;
                    
                    %w.r.t h
                    X_h=res;
                    for i =1:length(obj.stem_data.stem_crossval.stem_varset)
                        X_h(i,:,:)=obj.stem_data.stem_crossval.stem_varset{i}.X_h{1};
                    end
                    
                    tmp1=res(:);
                    tmp2=X_h(:);
                    tmp3=Y(:);
                    %
                    X_h_unique=unique(tmp2);
                    if length(X_h_unique)/length(tmp2)>0.3 
                        bins=linspace(min(tmp2),max(tmp2),21);
                        X_h_sel=[];
                        cut=sqrt(size(st_krig_data.grid.coordinate,1)*obj.stem_data.stem_datestamp.T);
                        for i=1:(length(bins)-1)
                            bin_h=tmp2(tmp2>=bins(i)&tmp2<bins(i+1));
                            if length(bin_h)>cut
                                X_h_sel=cat(1,X_h_sel,mean(bin_h));
                            end
                        end
                        if isempty(X_h_sel)
                            warning('Cross-validation MSE with respect to h can not be available due to sparse h.')
                        else
                            cv_mse_h=X_h_sel;
                            Y_h_var=X_h_sel;
                            for i=1:length(X_h_sel)
                                cv_mse_h(i) = nanvar(tmp1(tmp2>=bins(i)&tmp2<bins(i+1)));
                                Y_h_var(i)= nanvar(tmp3(tmp2>=bins(i)&tmp2<bins(i+1)));
                            end

                            obj.stem_crossval_result.cv_mse_h=cv_mse_h;
                            obj.stem_crossval_result.cv_R2_h=1-cv_mse_h./Y_h_var;
                            obj.stem_crossval_result.cv_h=X_h_sel;
                        end
                    else 
                        X_h_unique_sel=[];
                        for i=1:length(X_h_unique)
                            if sum(tmp2==X_h_unique(i))>sqrt(size(st_krig_data.grid.coordinate,1)*obj.stem_data.stem_datestamp.T)
                                X_h_unique_sel=cat(1,X_h_unique_sel,X_h_unique(i));
                            end
                        end
                        if isempty(X_h_unique_sel)
                            warning('Cross-validation MSE with respect to h can not be available due to sparse h.')
                        else
                            cv_mse_h=X_h_unique_sel;
                            Y_h_var=X_h_unique_sel;
                            for i=1:length(X_h_unique_sel)
                                cv_mse_h(i) = nanvar(tmp1(tmp2==X_h_unique_sel(i)));
                                Y_h_var(i)= nanvar(tmp3(tmp2==X_h_unique_sel(i)));
                            end

                            obj.stem_crossval_result.cv_mse_h=cv_mse_h;
                            obj.stem_crossval_result.cv_R2_h=1-cv_mse_h./Y_h_var;
                            obj.stem_crossval_result.cv_h=X_h_unique_sel;
                        end
                    end
                else
                    for i=1:length(st_krig_result)
                        if obj.stem_data.stem_crossval.stem_varset{i}.N>0
                            if iscell(st_krig_result)
                                obj.stem_crossval_result{i}.res=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}-st_krig_result{i}.y_hat;
                            else
                                obj.stem_crossval_result{i}.res=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}-st_krig_result.y_hat;
                            end

                            obj.stem_crossval_result{i}.cv_mse_s=(nanvar(obj.stem_crossval_result{i}.res'))';
                            obj.stem_crossval_result{i}.cv_mse_t=nanvar(obj.stem_crossval_result{i}.res);

                            obj.stem_crossval_result{i}.cv_R2_s = 1 - (obj.stem_crossval_result{i}.cv_mse_s'./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}'))';
                            obj.stem_crossval_result{i}.cv_R2_t = 1 - obj.stem_crossval_result{i}.cv_mse_t./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1});

                            %obj.stem_crossval_result{i}.mse=(nanvar(obj.stem_crossval_result{i}.res'))';
                            %obj.stem_crossval_result{i}.mse_time=nanvar(obj.stem_crossval_result{i}.res);

                            %obj.stem_crossval_result{i}.relative_mse=(obj.stem_crossval_result{i}.mse'./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}'))';
                            %obj.stem_crossval_result{i}.relative_mse_time=obj.stem_crossval_result{i}.mse_time./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1});

                            %obj.stem_crossval_result{i}.R2=1-obj.stem_crossval_result{i}.relative_mse;

                            if obj.stem_data.stem_varset_p.standardized
                                s=obj.stem_data.stem_varset_p.Y_stds{i};
                                m=obj.stem_data.stem_varset_p.Y_means{i};
                            else
                                s=1;
                                m=0;
                            end

                            if not(obj.stem_data.stem_varset_p.log_transformed)&&not(obj.stem_data.stem_varset_p.boxcox_transformed)
                                y_hat_back=st_krig_result{i}.y_hat*s+m;
                                y=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m;
                            end

                            if obj.stem_data.stem_varset_p.log_transformed
                                y_hat_back=st_krig_result{i}.y_hat;
                                var_y_hat=st_krig_result{i}.diag_Var_y_hat;
                                y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                                y=exp(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m);
                            end

                            if obj.stem_data.stem_varset_p.boxcox_transformed
                                lambda=obj.stem_data.stem_varset_p.Y_lambda{idx_var};

                                y_hat_back=st_krig_result{i}.y_hat;
                                var_y_hat=st_krig_result{i}.diag_Var_y_hat;

                                y_hat_back=(max(0,(lambda*(y_hat_back*s+m)+1)).^(1/lambda)).*(1+(var_y_hat*s^2*(1-lambda))./(2*(lambda*(y_hat_back*s+m)+1).^2));
                                y=(lambda*(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m)+1).^(1/lambda);
                            end

                            obj.stem_crossval_result{i}.res_back=y-y_hat_back;
                            obj.stem_crossval_result{i}.y_back=y;
                            obj.stem_crossval_result{i}.y_hat_back=y_hat_back;
                        end
                    end
                end
            else
                disp('The stem_model object does not include cross validation information');
            end
        end
        
        function data_summary(obj)
            %DESCRIPTION: print the information on the dataset
            %
            %INPUT
            %obj  - [stem_model object]   (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window  
            
            stem_misc.disp_star('D-STEM version')
            disp(' ');
            disp(['Model estimation is performed using ', stem_misc.ver]);
             
            stem_misc.disp_star('Brief data description')
            
            disp(['Number of variables: q=',num2str(obj.nvar)]);
            disp(' ');
            if obj.stem_data.stem_modeltype.is({'f-HDGM','HDGM'})
                disp('Point variables:');
                disp(['sum of points for each variable: Nq=',num2str(obj.stem_data.Np)]);
                disp(['sum of sites for each variable: N=',num2str(obj.stem_data.N/obj.nvar)])
                disp(' ');   
            else
                disp('Point variables:');
                if not(isempty(obj.stem_data.stem_varset_p))
                    disp(['Np (sum of points for each variable): ',num2str(obj.stem_data.Np)]);
                end
                disp(' '); 
                if not(isempty(obj.stem_data.stem_varset_b))
                    disp('Pixel variables:');
                    disp(['Nb (sum of pixels for each variable): ',num2str(obj.stem_data.Nb)]);
                    disp(' ');
                end
                disp(['N (sum of sites for each variable): ',num2str(obj.stem_data.N)])
                disp(' ');   
            end
            
            
            if isdatetime(obj.stem_data.stem_datestamp.date_start)
                disp('Date and time steps:')
                disp(['  Stating date  : ',datestr(obj.stem_data.stem_datestamp.date_start)]);
                disp(['  Ending date   : ',datestr(obj.stem_data.stem_datestamp.date_end)]);
            else
                disp('Time steps:')
            end
            disp(['  Temporal steps: ',num2str(obj.T)]);
            disp(' ');   
            
            disp('Bounding box of the point data:');
            if strcmp(obj.stem_data.stem_gridlist_p.grid{1}.unit,'deg')
                prefix_x='  Longitude ';
                prefix_y='  Latitude ';
                postfix=[' ',obj.stem_data.stem_gridlist_p.grid{1}.unit];
            else
                prefix_x='X ';
                prefix_y='Y ';
                postfix=[' ',obj.stem_data.stem_gridlist_p.grid{1}.unit];
            end
            disp([prefix_y,'min : ',num2str(obj.stem_data.stem_gridlist_p.box(1),'%05.2f'),postfix])
            disp([prefix_y,'max : ',num2str(obj.stem_data.stem_gridlist_p.box(2),'%05.2f'),postfix])
            disp([prefix_x,'min: ',num2str(obj.stem_data.stem_gridlist_p.box(3),'%05.2f'),postfix])
            disp([prefix_x,'max: ',num2str(obj.stem_data.stem_gridlist_p.box(4),'%05.2f'),postfix])
            disp(' ');

            if not(isempty(obj.stem_data.stem_varset_p.X_z_name))
                 disp(['Number of latent variables: p=', num2str(obj.stem_par.p)])
                 disp(' ')
            end
           
            if not(isempty(obj.stem_data.stem_varset_p.X_beta_name))
                disp('Covariates related to beta:');
                output=obj.stem_data.stem_varset_p.X_beta_name{1}';
                disp(output);
            end
          
        end
          
        function print(obj)
            %DESCRIPTION: print the information on the estimation result
            %
            %INPUT
            %obj  - [stem_model object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window          
            
            stem_misc.disp_star('Model estimation results')   
            if obj.estimated
                if obj.tapering
                    disp('* Tapering is enabled.');    
                    if not(isempty(obj.tapering_p))
                        disp(['  Point data tapering: ',num2str(obj.stem_data.stem_gridlist_p.tap),' km']);
                    else
                        disp('   Tapering is NOT enabled on point data');
                    end
                    if not(isempty(obj.tapering_b))
                        disp(['  Pixel data tapering: ',num2str(obj.stem_data.stem_gridlist_b.tap),' km']);
                    else
                        disp('  Tapering is NOT enabled on pixel data');
                    end
                else
                    disp('* Tapering is not enabled');
                end
                disp(' ');
                if not(isempty(obj.stem_EM_result.logL))
                    disp(['* Observed data log-likelihood: ',num2str(obj.stem_EM_result.logL,'%05.3f')]);
                else
                    disp('* Observed data log-likelihood: not computed. Use the method set_logL of the class stem_model.');
                end
                disp(' ');
                counter=1;
                counter_varcov=1;
                if not(isempty(obj.stem_par.beta))
                    if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            disp(['* Beta coefficients related to the point variable ',obj.stem_data.stem_varset_p.Y_name{i}]);
                            output=cell(length(obj.stem_data.stem_varset_p.X_beta_name{i})+1,3);
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            output{1,4}='|t|';
                            for j=1:length(obj.stem_data.stem_varset_p.X_beta_name{i})
                                output{j+1,1}=obj.stem_data.stem_varset_p.X_beta_name{i}{j};
                                output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                else
                                    output{j+1,3}='Not computed';
                                end
                                
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,4}=num2str(abs(obj.stem_par.beta(counter)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))),'%05.3f');
                                else
                                    output{j+1,4}='Not computed';
                                end
                                
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                            disp(output);
                        end
                        if not(isempty(obj.stem_data.stem_varset_b))
                            if not(isempty(obj.stem_data.stem_varset_b.X_beta))
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    disp(['* Beta coefficients related to the pixel variable ',obj.stem_data.stem_varset_b.Y_name{i}]);
                                    output=cell(length(obj.stem_data.stem_varset_b.X_beta_name{i})+1,3);
                                    output{1,1}='Loading coefficient';
                                    output{1,2}='Value';
                                    output{1,3}='Std';
                                    for j=1:length(obj.stem_data.stem_varset_b.X_beta_name{i})
                                        output{j+1,1}=obj.stem_data.stem_varset_b.X_beta_name{i}{j};
                                        output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                        if not(isempty(obj.stem_par.varcov))
                                            output{j+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                        else
                                            output{j+1,3}='Not computed';
                                        end
                                        counter=counter+1;
                                        counter_varcov=counter_varcov+1;
                                    end
                                    disp(output);
                                end
                            end
                        end
                    else
                        if obj.stem_par.flag_beta_spline ==1 
                            disp('* Beta Spline coefficients');
                            output=cell(length(obj.stem_data.X_beta_name{1})+1,4);
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            output{1,4}='|t|';
                            for j=1:length(obj.stem_data.X_beta_name{1})
                                output{j+1,1}=obj.stem_data.X_beta_name{1}{j};
                                output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                else
                                    output{j+1,3}='Not computed';
                                end
                                
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,4}=num2str(abs(obj.stem_par.beta(counter)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))),'%05.3f');
                                else
                                    output{j+1,4}='Not computed';
                                end
                                
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                            disp(output);
                        else  
                            disp('* Beta coefficients');
                            output=cell(length(obj.stem_data.stem_varset_p.X_beta_name{1})+1,4);
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            output{1,4}='|t|';
                            for j=1:length(obj.stem_data.stem_varset_p.X_beta_name{1})
                                output{j+1,1}=obj.stem_data.stem_varset_p.X_beta_name{1}{j};
                                output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                else
                                    output{j+1,3}='Not computed';
                                end
                                
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,4}=num2str(abs(obj.stem_par.beta(counter)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))),'%05.3f');
                                else
                                    output{j+1,4}='Not computed';
                                end
                                
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                            disp(output);
                        end   
                    end
                end
                
                if obj.stem_par.stem_par_constraints.sigma_eps_fixed==0
                    if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                        output=cell(obj.stem_data.stem_varset_p.nvar+1,3);
                        if size(obj.stem_par.sigma_eps,3)==1
                            disp('* Sigma^2_eps')
                            output{1,1}='Variable';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            for i=1:obj.stem_data.stem_varset_p.nvar
                                output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                                output{i+1,2}=num2str(obj.stem_par.sigma_eps(i,i),'%05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                else
                                    output{i+1,3}='Not computed';
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                            if not(isempty(obj.stem_data.stem_varset_b))
                                delta=obj.stem_data.stem_varset_p.nvar;
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                                    output{i+1+delta,2}=num2str(obj.stem_par.sigma_eps(i+delta,i+delta),'%05.3f');
                                    if not(isempty(obj.stem_par.varcov))
                                        output{i+1+delta,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                    else
                                        output{i+1+delta,3}='Not computed';
                                    end
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                        else
                            disp('* Sigma^2_eps (time average)')
                            T=size(obj.stem_par.sigma_eps,3);
                            step=obj.stem_data.stem_varset_p.nvar;
                            if not(isempty(obj.stem_data.stem_varset_b))
                                step=step+obj.stem_data.stem_varset_b.nvar;
                            end
                            output{1,1}='Variable';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            counter_varcov_temp=counter_varcov;
                            for i=1:obj.stem_data.stem_varset_p.nvar
                                output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                                output{i+1,2}=num2str(mean(obj.stem_par.sigma_eps(i,i,:)),'%05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,3}=num2str(mean(sqrt(diag(obj.stem_par.varcov(counter_varcov_temp+(i-1):step:counter_varcov_temp+step*T-1,counter_varcov_temp+(i-1):step:counter_varcov_temp+step*T-1)))),'%05.3f');
                                else
                                    output{i+1,3}='Not computed';
                                end
                                counter=counter+1;
                            end
                            counter_varcov=counter_varcov+T*obj.stem_data.stem_varset_p.nvar;
                            if not(isempty(obj.stem_data.stem_varset_b))
                                delta=obj.stem_data.stem_varset_p.nvar;
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                                    output{i+1+delta,2}=num2str(mean(obj.stem_par.sigma_eps(i+delta,i+delta,:)),'%05.3f');
                                    if not(isempty(obj.stem_par.varcov))
                                        output{i+1+delta,3}=num2str(sqrt(diag(obj.stem_par.varcov(counter_varcov_temp+(i+delta-1):step:counter_varcov_temp+step*T-1,counter_varcov_temp+(i+delta-1):step:counter_varcov_temp+step*T-1))),'%05.3f');
                                    else
                                        output{i+1+delta,3}='Not computed';
                                    end
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                                counter_varcov=counter_varcov+T*obj.stem_data.stem_varset_b.nvar;
                            end
                        end
                        
                        disp(output);
                    else
                        if obj.stem_par.flag_sigma_eps_spline==1 
                            output=cell(1+obj.stem_par.k_sigma,4);
                            disp('* log(Sigma_eps) spline coefficients')
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            output{1,4}='|t|';
                            for j=1:obj.stem_par.k_sigma
                                output{j+1,1}=['Basis_',num2str(j)];
                                output{j+1,2}=num2str(obj.stem_par.sigma_eps(j),'%+05.3f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                                else
                                    output{j+1,3}='Not computed';
                                end
                                if not(isempty(obj.stem_par.varcov))
                                    output{j+1,4}=num2str(abs(obj.stem_par.sigma_eps(j)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))),'%05.3f');
                                else
                                    output{j+1,4}='Not computed';
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end 
                            disp(output);
                        else
                            output=cell(2,2);
                            disp('* Sigma_eps diagonal elements (Variance)')
                            output{1,1}='Value';
                            output{1,2}='Std';
                            output{2,1}=num2str(obj.stem_par.sigma_eps,'%05.3f');
                            if not(isempty(obj.stem_par.varcov))
                                output{2,2}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                            else
                                output{2,2}='Not computed';
                            end
                            counter=counter+1;
                            counter_varcov=counter_varcov+1;
                            disp(output);
                        end
                    end
                end
                if not(isempty(obj.stem_data.stem_varset_b))
                    disp('* alpha_bp elements')
                    output=cell(obj.stem_data.stem_varset_p.nvar+1,3);
                    output{1,1}='Variable';
                    output{1,2}='Value';
                    output{1,3}='Std';
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        output{i+1,2}=num2str(obj.stem_par.alpha_bp(i),'%+05.3f');
                        if not(isempty(obj.stem_par.varcov))
                            output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                        else
                            output{i+1,3}='Not computed';
                        end
                        counter=counter+1;
                        counter_varcov=counter_varcov+1;
                    end
                   
                    disp(output);
                 
                    output=[];
                    if obj.stem_par.stem_par_constraints.pixel_correlated
                        disp('* Pixel data are cross-correlated.');
                        disp(' ');
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        output{2,1}='Theta_b';
                        output{2,2}=num2str(obj.stem_par.theta_b(1),'%05.3f');
                        if not(isempty(obj.stem_par.varcov))
                            output{2,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                        else
                            output{2,3}='Not computed';
                        end
                        counter=counter+1;
                        counter_varcov=counter_varcov+1;
                        disp(output);
                        output=cell(obj.stem_data.stem_varset_b.nvar+1,obj.stem_data.stem_varset_b.nvar+1);
                        disp('* v_b matrix:');
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            for j=i+1:obj.stem_data.stem_varset_b.nvar
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_b(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_b(i,j),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                        end
                        disp(output);
                    else
                        disp('* Pixel data are NOT cross-correlated.');
                        disp(' ');
                        disp('* Theta_b elements:');
                        output=cell(obj.stem_data.stem_varset_b.nvar+1,3);
                        output{1,1}='Variable';
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,2}=num2str(obj.stem_par.theta_b(i),'%05.3f');
                            if not(isempty(obj.stem_par.varcov))
                                output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f');
                            else
                                output{i+1,3}='Not computed';
                            end
                            counter=counter+1;
                            counter_varcov=counter_varcov+1;
                        end
                        disp(output);
                    end
                end
                output=cell(obj.stem_data.stem_varset_p.nvar*2,obj.stem_par.k+1);
                if obj.stem_par.k>0
                    disp(['* ',num2str(obj.stem_par.k),' fine-scale coregionalization components w_p']);
                    disp(' ');
                    disp('* theta_p elements:');
                    output=cell(1,3);
                    output{1,1}='Coreg. component';
                    output{1,2}='Value [km]';
                    output{1,3}='Std [km]';
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end
                        output{k+1,1}=[num2str(k),postfix];
                        output{k+1,2}=num2str(obj.stem_par.theta_p(:,k),'%06.2f');
                        if not(isempty(obj.stem_par.varcov))
                            output{k+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f');
                        else
                            output{k+1,3}='Not computed';
                        end
                        counter=counter+1;
                        counter_varcov=counter_varcov+1;
                    end
                    disp(output);
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end                        
                        disp(['* v_p matrix for the ',num2str(k),postfix,' coreg. component:']);
                        output=cell(obj.stem_data.stem_varset_p.nvar+1,obj.stem_data.stem_varset_p.nvar+1);
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        end
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=i:obj.stem_data.stem_varset_p.nvar
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3))),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_p(i,j,min(k,size(obj.stem_par.v_p,3))),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                        end                        
                        disp(output);
                    end
                end
                
                if obj.stem_par.p>0
                    output=cell(obj.stem_par.p+1,obj.stem_par.p+1);
                    disp('* Transition matrix G:');
                    c=1;
                    if obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                        if obj.stem_data.stem_modeltype.is({'f-HDGM'})
                            output=cell(obj.stem_par.p+1,3);
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            for j=1:obj.stem_par.p
                                output{c+1,1}=obj.stem_data.stem_varset_p.X_z_name{1}{j};
                                c=c+1;
                            end
                        elseif isdiag(obj.stem_par.G)
                            output=cell(obj.stem_par.p+1,3);
                            output{1,1}='Loading coefficient';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            for j=1:obj.stem_par.p
                                output{c+1,1}=obj.stem_data.stem_varset_p.X_z_name{j}{:};
                                c=c+1;
                            end
                        else
                            for j=1:obj.stem_par.p
                                output{1,c+1}=obj.stem_data.stem_varset_p.X_z_name{1}{j};
                                output{c+1,1}=output{1,c+1};
                                c=c+1;
                            end
                        end  
                    else
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                                output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                                output{c+1,1}=output{1,c+1};
                                c=c+1;
                            end
                        end
                        if not(isempty(obj.stem_data.stem_varset_b))
                            if not(isempty(obj.stem_data.stem_varset_b.X_z))
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                        output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                        output{c+1,1}=output{1,c+1};
                                        c=c+1;
                                    end
                                end
                            end
                        end
                    end
                    
                    if obj.stem_par.stem_par_constraints.time_diagonal || obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                        if isdiag(obj.stem_par.G) && obj.stem_data.stem_modeltype.is({'f-HDGM'})
                            temp=abs(obj.stem_par.G-eye(size(obj.stem_par.G)));
                            if sum(temp(:))==0
                                random_walk=1;
                            else
                                random_walk=0;
                            end
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_par.varcov))&&random_walk==0
                                    output{i+1,2}=num2str(obj.stem_par.G(i,i),'%+05.3f');
                                    output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.3f'); 
                                else
                                    output{i+1,2}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                end
                                counter=counter+1;
                                if not(random_walk)
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                        else
                            temp=abs(obj.stem_par.G-eye(size(obj.stem_par.G)));
                            if sum(temp(:))==0
                                random_walk=1;
                            else
                                random_walk=0;
                            end
                            for j=1:obj.stem_par.p
                                for i=1:obj.stem_par.p
                                    if i==j
                                        if not(isempty(obj.stem_par.varcov))&&random_walk==0
                                            if isdiag(obj.stem_par.G)
                                                output{i+1,2}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                                output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter,counter)),'%03.2f');
                                            else
                                                output{i+1,i+1}=[num2str(obj.stem_par.G(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                            end 
                                        else
                                            output{i+1,i+1}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                        end
                                        counter=counter+1;
                                        if not(random_walk)
                                            counter_varcov=counter_varcov+1;
                                        end
                                    else
                                        output{i+1,j+1}='0';
                                    end
                                end
                            end
                        end  
                    else
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.G(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.G(i,j),'%+05.2f');
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                        end
                    end
                    disp(output);
                    
                    if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                        disp('* Sigma_eta matrix:');
                        c=1;
                        output=cell(obj.stem_par.p+1,obj.stem_par.p);
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                                output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                                output{c+1,1}=output{1,c+1};
                                c=c+1;
                            end
                        end
                        if not(isempty(obj.stem_data.stem_varset_b))
                            if not(isempty(obj.stem_data.stem_varset_b.X_z))
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                        output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                        output{c+1,1}=output{1,c+1};
                                        c=c+1;
                                    end
                                end
                            end
                        end
                        if obj.stem_par.stem_par_constraints.time_diagonal
                            for j=1:obj.stem_par.p
                                for i=1:obj.stem_par.p
                                    if i==j
                                        if not(isempty(obj.stem_par.varcov))
                                            output{i+1,i+1}=[num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                        else
                                            output{i+1,i+1}=num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f');
                                        end
                                        counter=counter+1;
                                        counter_varcov=counter_varcov+1;
                                    else
                                        output{i+1,j+1}='0';
                                    end
                                end
                            end
                        else
                            for j=1:obj.stem_par.p
                                for i=j:obj.stem_par.p
                                    if not(isempty(obj.stem_par.varcov))
                                        output{i+1,j+1}=[num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                    else
                                        output{i+1,j+1}=num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f');
                                    end
                                    output{j+1,i+1}=output{i+1,j+1};
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                        end
                        disp(output);
                    else
                        disp('* Fine-scale coregionalization components z');
                        disp(' ');
                        
                        disp('* theta_z elements:');
                        if not(strcmp(obj.stem_par.correlation_type,'expsphere'))
                            if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                                output=cell(length(obj.stem_par.theta_z)+1,2);
                                output{1,1}='Value [km]';
                                output{1,2}='Std [km]';
                                output{2,1}=num2str(obj.stem_par.theta_z,'%06.2f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{2,2}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f');
                                else
                                    output{2,2}='Not computed';
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            else
                                output=cell(length(obj.stem_par.theta_z)+1,3);
                                output{1,1}='Basis';
                                output{1,2}='Value [km]';
                                output{1,3}='Std [km]';
                                for b=1:length(obj.stem_par.theta_z)
                                    output{b+1,1}=obj.stem_data.stem_varset_p.X_z_name{1}{b};
                                    output{b+1,2}=num2str(obj.stem_par.theta_z(b),'%06.2f');
                                    if not(isempty(obj.stem_par.varcov))
                                        output{b+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f');
                                    else
                                        output{b+1,3}='Not computed';
                                    end
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                        else
                            if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                                output=cell(3,3);
                                output{1,1}='Name';
                                output{1,2}='Value';
                                output{1,3}='Std';
                                output{2,1}='Spatial range (km)';
                                output{3,1}='Anisotropy range (deg)';
                                
                                output{2,2}=num2str(obj.stem_par.theta_z(1),'%06.2f');
                                output{3,2}=num2str(obj.stem_par.theta_z(2),'%06.2f');
                                if not(isempty(obj.stem_par.varcov))
                                    output{2,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f');
                                    output{3,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov+1,counter_varcov+1)),'%05.2f');
                                else
                                    output{2,3}='Not computed';
                                    output{3,3}='Not computed';
                                end
                                counter=counter+2;
                                counter_varcov=counter_varcov+2;
                            else
                                output=cell(1+length(obj.stem_par.theta_z),5);
                                output{1,1}='Basis';
                                output{1,2}='Spat. range (km)';
                                output{1,3}='Anis. range (km)';
                                output{1,4}='Spat. range std (km)';
                                output{1,5}='Anis. range std (km)';
                                for b=1:length(obj.stem_par.theta_z)
                                    output{b+1,1}=obj.stem_data.stem_varset_p.X_z_name{1, 1}{1, b};
                                    output{b+1,2}=num2str(obj.stem_par.theta_z(1,b),'%06.2f');
                                    output{b+1,3}=num2str(obj.stem_par.theta_z(2,b),'%06.2f');
                                    if not(isempty(obj.stem_par.varcov))
                                        output{b+1,4}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%05.2f');
                                        output{b+1,5}=num2str(sqrt(obj.stem_par.varcov(counter_varcov+1,counter_varcov+1)),'%05.2f');
                                    else
                                        output{b+1,4}='Not computed';
                                        output{b+1,5}='Not computed';
                                    end
                                    counter=counter+2;
                                    counter_varcov=counter_varcov+2;
                                end
                            end
                        end
                        disp(output);
                        
                        if obj.stem_data.stem_modeltype.is('HDGM')
                            disp('* v_z matrix for the coreg. component:');
                            output=cell(obj.stem_par.p+1,obj.stem_par.p+1);
                            for i=1:obj.stem_par.p
                                output{1,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                                output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                            end
                            for i=1:obj.stem_par.p
                                for j=i:obj.stem_par.p
                                    if not(isempty(obj.stem_par.varcov))
                                        output{i+1,j+1}=[num2str(obj.stem_par.v_z(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f'),')'];
                                        output{j+1,i+1}=output{i+1,j+1};
                                    else
                                        output{i+1,j+1}=num2str(obj.stem_par.v_z(i,j),'%+05.2f');
                                        output{j+1,i+1}=output{i+1,j+1};
                                    end
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                            disp(output);
                        end
                        if obj.stem_data.stem_modeltype.is('f-HDGM')
                            disp('* v_z matrix for the coreg. component:');
                            
                            output=cell(obj.stem_par.p+1,3);
                            output{1,1}='Basis';
                            output{1,2}='Value';
                            output{1,3}='Std';
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_par.varcov))
                                    output{i+1,1}=obj.stem_data.stem_varset_p.X_z_name{1, 1}{1, i};
                                    output{i+1,2}=num2str(obj.stem_par.v_z(i,i),'%+05.2f');
                                    output{i+1,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)),'%03.2f');
                                else
                                    output{i+1,1}=obj.stem_data.stem_varset_p.X_z_name{1, 1}{1, i};
                                    output{i+1,2}=num2str(obj.stem_par.v_z(i,i),'%+05.2f');
                                    output{i+1,3}='Not computed';
                                end
                                counter=counter+1;
                                counter_varcov=counter_varcov+1;
                            end
                            disp(output);
                        end
                    end
                end
                if obj.stem_data.stem_modeltype.is('f-HDGM')
                    if strcmp(obj.stem_data.stem_fda.spline_type, 'Fourier')
                        disp(['* The basis number for z is ', num2str(obj.stem_data.stem_fda.spline_nbasis_z)])
                        if obj.stem_par.flag_beta_spline ==1
                            disp(['* The basis number for beta(h) is ', num2str(obj.stem_data.stem_fda.spline_nbasis_beta)])
                        end
                        if obj.stem_par.flag_sigma_eps_spline ==1
                            disp(['* The basis number for sigma_eps^2(h) is ', num2str(obj.stem_data.stem_fda.spline_nbasis_sigma)])
                        end
                    elseif strcmp(obj.stem_data.stem_fda.spline_type, 'Bspline')
                        disp(['* The knots for z is ', num2str(obj.stem_data.stem_fda.spline_knots_z)])
                        if obj.stem_par.flag_beta_spline ==1
                            disp(['* The knots for beta(h£© is ', num2str(obj.stem_data.stem_fda.spline_knots_beta)])
                        end
                        if obj.stem_par.flag_sigma_eps_spline ==1
                            disp(['* The knots for sigma_eps^2(h) is ',num2str(obj.stem_data.stem_fda.spline_knots_sigma)])
                        end
                    end
                end
            else
                disp('The model has not been estimated yet. Use the method print of the class stem_data to print data information.');
            end            
        end
        
        function plot_par(obj,vertical)
            %DESCRIPTION: plot the beta(h) and sigma_eps2(h) when modeltype is f-HDGM
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %vertical - [boolean]          (1x1) indicating if the plot is vertical or not.
            %
            %OUTPUT
            %
            %None: the plot of beta(h) and sigma_eps2(h) is in the figure window
            
            if nargin==1
                vertical=0;
            end
            
            if not(obj.stem_par.stem_modeltype.is('f-HDGM'))
                error('The parplot function is only used when model type is f-HDGM')
            end
            if (~obj.stem_par.flag_beta_spline)&&(~obj.stem_par.flag_sigma_eps_spline)
                warning('Since the number of basis to estimate beta(h) and sigma_eps2(h) is none, nothing is to plot')
            end
            if isempty(obj.stem_par.varcov)
                warning('Variance covariance matrix is not estimated')
            end
            Num = 0;
            if obj.stem_par.flag_beta_spline==1
                Num = Num + length(obj.stem_data.stem_varset_p.X_beta_name{1,1});
            end
            if obj.stem_par.flag_sigma_eps_spline==1
                Num = Num + 1;
            end
            nrow = floor(sqrt(Num));
            ncol = ceil(Num/nrow);
            
            quant = [0.005, 0.025, 0.05];
            coef = -norminv(quant,0,1);

            counter = 0;
            
            figure
            if vertical
                if obj.stem_par.flag_beta_spline==1
                    disp('****************');
                    disp('\Beta_(h) plot');
                    disp('****************');
                    range = obj.stem_data.stem_fda.spline_range;
                    h = range(1):0.1:range(2);
                    k = obj.stem_par.k_beta;
                    if not(isempty(obj.stem_par.varcov))
                        for i = 1:length(obj.stem_data.stem_varset_p.X_beta_name{1,1})
                            covariate = ['\beta_{',obj.stem_data.stem_varset_p.X_beta_name{1,1}{i}, '}(',obj.stem_data.stem_varset_p.X_h_name,')'];
                            basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_beta));
                            beta = obj.stem_par.beta((i-1)*k+(1:k));
                            beta_h = basis*beta;
                            subplot(nrow,ncol,i);
                            plot(beta_h, h,'k');
                            for j=1:3
                                beta_h_up = beta_h + coef(j)*sqrt(diag(basis*obj.stem_par.varcov((i-1)*k+(1:k),(i-1)*k+(1:k))*basis'));
                                beta_h_low = beta_h - coef(j)*sqrt(diag(basis*obj.stem_par.varcov((i-1)*k+(1:k),(i-1)*k+(1:k))*basis'));

                                patch([beta_h_up' fliplr(beta_h_low')],[h fliplr(h)],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                            end
                            title(covariate,'FontSize',14);
                            set(gca,'Ydir','reverse','FontSize',14); 
                            ylim([h(1),h(end)]);
                            ylabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        end
                        counter = counter + length(obj.stem_data.X_beta_name{1,1});
                    else
                        for i = 1:length(obj.stem_data.stem_varset_p.X_beta_name{1,1})
                            covariate = ['\beta_{',obj.stem_data.stem_varset_p.X_beta_name{1,1}{i}, '}(',obj.stem_data.stem_varset_p.X_h_name,')'];
                            basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_beta));
                            beta = obj.stem_par.beta((i-1)*k+(1:k));
                            beta_h = basis*beta;
                            subplot(nrow,ncol,i);
                            plot(beta_h,h,'k');
                            title(covariate,'FontSize',14);
                            set(gca,'Ydir','reverse','FontSize',14); 
                            xlabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        end
                        counter = counter + length(obj.stem_data.X_beta_name{1,1});
                    end 
                else
                    counter=length(obj.stem_data.stem_varset_p.X_beta_name{1,1});
                end

                if obj.stem_par.flag_sigma_eps_spline==1
                    disp('****************');
                    disp('Sigma(h) plot');
                    disp('****************');
                    range = obj.stem_data.stem_fda.spline_range;
                    h = range(1):0.1:range(2);
                    k = obj.stem_par.k_sigma;
                    if not(isempty(obj.stem_par.varcov))
                        covariate = ['\sigma_{\epsilon}^2(',obj.stem_data.stem_varset_p.X_h_name,')'];
                        basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_sigma));
                        sigma_eps = obj.stem_par.sigma_eps;
                        if obj.stem_par.flag_logsigma==1
                            sigma_eps_h = exp(basis*sigma_eps);
                            Varcov_eps = obj.stem_par.varcov(counter+(1:k),counter+(1:k));
                            Varcov = diag(basis*Varcov_eps*basis');
                            Varcov1 = (exp(Varcov)-ones(length(Varcov),1)).*((sigma_eps_h).^2).*(exp(Varcov));

                            if obj.stem_par.flag_beta_spline==1
                                subplot(nrow,ncol,length(obj.stem_data.stem_varset_p.X_beta_name{1,1})+1);
                            else
                                subplot(nrow,ncol,1);
                            end
                            plot(sigma_eps_h,h,'k');
                            for j=1:3
                                sigma_eps_h_up = sigma_eps_h + coef(j)*sqrt(Varcov1);
                                sigma_eps_h_low = sigma_eps_h - coef(j)*sqrt(Varcov1); 
                                patch([sigma_eps_h_up' fliplr(sigma_eps_h_low')],[h fliplr(h)],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                            end
                            ylim([h(1), h(end)]);
                            title(covariate,'FontSize',14);
                            set(gca,'Ydir','reverse','FontSize',14);
                            ylabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        end
                    else
                        covariate = ['\sigma_{\epsilon}^2(',obj.stem_data.stem_varset_p.X_h_name,')'];
                        basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_sigma));
                        sigma_eps = obj.stem_par.sigma_eps;
                        if obj.stem_par.flag_logsigma==1
                            sigma_eps_h = exp(basis*sigma_eps);
                            if obj.stem_par.flag_beta_spline==1
                                subplot(nrow,ncol,length(obj.stem_data.stem_varset_p.X_beta_name{1,1})+1);
                            else
                                subplot(nrow,ncol,1);
                            end
                            plot(sigma_eps_h,h,'k');
                            title(covariate,'FontSize',14);
                            set(gca,'Ydir','reverse','FontSize',14);
                            ylabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        end 
                    end
                end
            else
                if obj.stem_par.flag_beta_spline==1
                    disp('****************');
                    disp('\Beta_(f) plot');
                    disp('****************');
                    range = obj.stem_data.stem_fda.spline_range;
                    h = range(1):0.1:range(2);
                    k = obj.stem_par.k_beta;
                    if not(isempty(obj.stem_par.varcov))
                        for i = 1:length(obj.stem_data.stem_varset_p.X_beta_name{1,1})
                            covariate = ['\beta_{',obj.stem_data.stem_varset_p.X_beta_name{1,1}{i}, '}(',obj.stem_data.stem_varset_p.X_h_name,')'];
                            basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_beta));
                            beta = obj.stem_par.beta((i-1)*k+(1:k));
                            beta_h = basis*beta;
                            subplot(nrow,ncol,i);
                            plot(h, beta_h, 'k');
                            for j=1:3
                                beta_h_up = beta_h + coef(j)*sqrt(diag(basis*obj.stem_par.varcov((i-1)*k+(1:k),(i-1)*k+(1:k))*basis'));
                                beta_h_low = beta_h - coef(j)*sqrt(diag(basis*obj.stem_par.varcov((i-1)*k+(1:k),(i-1)*k+(1:k))*basis'));

                                patch([h fliplr(h)],[beta_h_up' fliplr(beta_h_low')],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                            end
                            title(covariate,'FontSize',14);
                            set(gca,'FontSize',14);
                            xlabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                            xlim([h(1),h(end)]);
                        end
                        counter = counter + length(obj.stem_data.X_beta_name{1,1});
                    else
                        for i = 1:length(obj.stem_data.stem_varset_p.X_beta_name{1,1})
                            covariate = ['\beta_{',obj.stem_data.stem_varset_p.X_beta_name{1,1}{i}, '}(',obj.stem_data.stem_varset_p.X_h_name,')'];
                            basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_beta));
                            beta = obj.stem_par.beta((i-1)*k+(1:k));
                            beta_h = basis*beta;
                            subplot(nrow,ncol,i);
                            plot(h, beta_h,'k');
                            title(covariate,'FontSize',14);
                            set(gca,'FontSize',14);
                            xlabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                            xlim([h(1),h(end)]);
                        end
                        counter = counter + length(obj.stem_data.X_beta_name{1,1});
                    end 
                else
                    counter=length(obj.stem_data.stem_varset_p.X_beta_name{1,1});
                end

                if obj.stem_par.flag_sigma_eps_spline==1
                    disp('****************');
                    disp('Sigma(h) plot');
                    disp('****************');
                    range = obj.stem_data.stem_fda.spline_range;
                    h = range(1):0.1:range(2);
                    k = obj.stem_par.k_sigma;
                    if not(isempty(obj.stem_par.varcov))
                        covariate = ['\sigma_{\epsilon}^2(',obj.stem_data.stem_varset_p.X_h_name,')'];
                        basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_sigma));
                        sigma_eps = obj.stem_par.sigma_eps;
                        if obj.stem_par.flag_logsigma==1
                            sigma_eps_h = exp(basis*sigma_eps);
                            Varcov_eps = obj.stem_par.varcov(counter+(1:k),counter+(1:k));
                            Varcov = diag(basis*Varcov_eps*basis');
                            Varcov1 = (exp(Varcov)-ones(length(Varcov),1)).*((sigma_eps_h).^2).*(exp(Varcov));
                            if obj.stem_par.flag_beta_spline==1
                                subplot(nrow,ncol,length(obj.stem_data.stem_varset_p.X_beta_name{1,1})+1);
                            else
                                subplot(nrow,ncol,1);
                            end
                            plot(h, sigma_eps_h,'k');
                            for j=1:3
                                sigma_eps_h_up = sigma_eps_h + coef(j)*sqrt(Varcov1);
                                sigma_eps_h_low = sigma_eps_h - coef(j)*sqrt(Varcov1);  
                                patch([h fliplr(h)],[sigma_eps_h_up' fliplr(sigma_eps_h_low')],'r','FaceAlpha',0.15,'EdgeAlpha',0)
                            end
                            xlim([h(1), h(end)]);
                            title(covariate,'FontSize',14);
                            set(gca,'FontSize',14);
                            xlabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        else 
                            %Currently sigma_eps(h) is estimated by log
                            %transmission to enfore the positive value of
                            %sigma_eps.
                            %Nothing to do here.
                        end 
                    else
                        covariate = ['\sigma_{\epsilon}^2(',obj.stem_data.stem_varset_p.X_h_name,')'];
                        basis = full(getbasismatrix(h,obj.stem_data.stem_fda.spline_basis_sigma));
                        sigma_eps = obj.stem_par.sigma_eps;
                        if obj.stem_par.flag_logsigma==1
                            sigma_eps_h = exp(basis*sigma_eps);
                            if obj.stem_par.flag_beta_spline==1
                                subplot(nrow,ncol,length(obj.stem_data.stem_varset_p.X_beta_name{1,1})+1);
                            else
                                subplot(nrow,ncol,1);
                            end
                            plot(h, sigma_eps_h,'k');
                            xlim([h(1), h(end)]);
                            title(covariate,'FontSize',14);
                            set(gca,'FontSize',14);
                            xlabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
                        else
                            %Nothing to do here.
                        end 
                    end
                end
            end
           
        end 
        
        function plot_xval(obj,vertical,variable_name)
            %DESCRIPTION: the cross-validation results when modeltype is f-HDGM
            %
            %INPUT
            %obj              - [stem_model object] (1x1) the stem_model object
            %vertical         - [boolean]           (1x1) indicating if the plot is vertical or not.
            %<variable_name>  -[string]              (1x1) the names of the cross-validation variables
            %
            %OUTPUT
            %
            %None: the cross-validation results are shown in the figure window
            
            if not(obj.stem_par.stem_modeltype.is('f-HDGM')||obj.stem_par.stem_modeltype.is('HDGM'))
                error('The plot_xval method is only used when model type is f-HDGM')
            end
            
            if not(obj.cross_validation)
                error('The plot_xval method is only used when crossvalidation is implemented')
            end
            
            if nargin==1
                vertical=0;
                variable_name=[];
            end
            
            if nargin==2
                variable_name=[];
            end
            
            if abs(obj.stem_data.stem_gridlist_p.box(4)-obj.stem_data.stem_gridlist_p.box(3))>180
                global_idx=1;
            else
                global_idx=0;
            end
            
            if isempty(variable_name)
                idx=1;
            else
                idx=obj.stem_data.stem_varset_p.get_Y_index(variable_name);
                if isempty(idx)
                    error('Variable not found');
                end
            end
                
            if iscell(obj.stem_crossval_result) 
                xvalresult=obj.stem_crossval_result{idx};
            else
                xvalresult=obj.stem_crossval_result;
            end
            
            
            % plot of mse_s and r2_s
            if global_idx
                figure 
                temp=xvalresult.cv_mse_s;
                lat=obj.stem_data.stem_crossval.stem_gridlist{idx}.grid{1}.coordinate(:,1);
                lon=obj.stem_data.stem_crossval.stem_gridlist{idx}.grid{1}.coordinate(:,2);

                scatter(lon,lat,100,temp,'filled');
                cl = colorbar;
                cl.Limits=[min(temp(:)) max(temp(:))];

                lat1=obj.stem_data.stem_gridlist_p.grid{idx}.coordinate(:,1);
                lon1=obj.stem_data.stem_gridlist_p.grid{idx}.coordinate(:,2);
                geoshow(lat1,lon1,'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
                
                if not(isempty(obj.stem_data.shape))
                    try
                        geoshow(obj.stem_data.shape,'FaceColor','none');
                    catch
                        geoshow(obj.stem_data.shape);
                    end
                end

                if strcmp(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'deg')
                    xlabel('Longitude [deg]','FontSize',14);
                    ylabel('Latitude [deg]','FontSize',14);
                else
                    xlabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                    ylabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                end 
                title('CV-MSE(s)','FontSize',14,'interpreter','latex');
                set(gca,'FontSize',14);
                set(gcf, 'renderer', 'zbuffer');
                xticks(-180:30:180);
                yticks(-90:15:90);
                xlim([-180, 180]);
                ylim([-90, 90]);
                box on
                grid on

                figure 
                temp=xvalresult.cv_R2_s;
                scatter(lon,lat,100,temp,'filled');

                geoshow(lat1,lon1,'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
               
                if not(isempty(obj.stem_data.shape))
                    try
                        geoshow(obj.stem_data.shape,'FaceColor','none');
                    catch
                        geoshow(obj.stem_data.shape);
                    end
                end

                if strcmp(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'deg')
                    xlabel('Longitude [deg]','FontSize',14);
                    ylabel('Latitude [deg]','FontSize',14);
                else
                    xlabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                    ylabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                end
                title('CV-R2(s)','FontSize',14,'interpreter','latex');
                cl = colorbar;
                cl.Limits=[min(temp(:)) max(temp(:))];
                set(gca,'FontSize',14);
                set(gcf, 'renderer', 'zbuffer');
                xticks(-180:30:180);
                yticks(-90:15:90);
                xlim([-180, 180]);
                ylim([-90, 90]);
                box on
                grid on 
            else
                figure          
                subplot(1,2,1)
                temp=xvalresult.cv_mse_s;
                lat=obj.stem_data.stem_crossval.stem_gridlist{idx}.grid{1}.coordinate(:,1);
                lon=obj.stem_data.stem_crossval.stem_gridlist{idx}.grid{1}.coordinate(:,2);

                scatter(lon,lat,100,temp,'filled');
                cl = colorbar;
                cl.Limits=[min(temp(:)) max(temp(:))];

                lat1=obj.stem_data.stem_gridlist_p.grid{idx}.coordinate(:,1);
                lon1=obj.stem_data.stem_gridlist_p.grid{idx}.coordinate(:,2);
                geoshow(lat1,lon1,'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
                
                if strcmp(obj.stem_data.shape(1,1).Geometry,'Line')
                    flag_readshape_from_user=1;
                else
                    flag_readshape_from_user=0;
                end
                
                if flag_readshape_from_user
                    if not(isempty(obj.stem_data.shape))
                        try
                            geoshow(obj.stem_data.shape,'FaceColor','none');
                        catch
                            geoshow(obj.stem_data.shape);
                        end
                    end
                end
                
                axis tight

                if strcmp(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'deg')
                    xlabel('Longitude [deg]','FontSize',14);
                    ylabel('Latitude [deg]','FontSize',14);
                else
                    xlabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                    ylabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                end 

                title('CV-MSE(s)','FontSize',14,'interpreter','latex');
                set(gca,'FontSize',14);
                box on
                grid on

                subplot(1,2,2)
                temp=xvalresult.cv_R2_s;
                scatter(lon,lat,100,temp,'filled');

                geoshow(lat1,lon1,'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
                
                if strcmp(obj.stem_data.shape(1,1).Geometry,'Line')
                    flag_readshape_from_user=1;
                else
                    flag_readshape_from_user=0;
                end
                
                if flag_readshape_from_user
                    if not(isempty(obj.stem_data.shape))
                        try
                            geoshow(obj.stem_data.shape,'FaceColor','none');
                        catch
                            geoshow(obj.stem_data.shape);
                        end
                    end
                end
                
                axis tight

                if strcmp(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'deg')
                    xlabel('Longitude [deg]','FontSize',14);
                    ylabel('Latitude [deg]','FontSize',14);
                else
                    xlabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                    ylabel(obj.stem_data.stem_gridlist_p.grid{1,1}.unit,'FontSize',14);
                end
                title('CV-R2(s)','FontSize',14,'interpreter','latex');
                cl = colorbar;
                cl.Limits=[min(temp(:)) max(temp(:))];
                set(gca,'FontSize',14);
                box on
                grid on  
            end
            
            if obj.stem_par.stem_modeltype.is('f-HDGM')
                % X_h_num for two cases
                X_h=nan(obj.stem_par.q,size(obj.stem_crossval_result.res{1},1),obj.stem_data.T);
                for i =1:length(obj.stem_data.stem_crossval.stem_varset)
                    X_h(i,:,:)=obj.stem_data.stem_crossval.stem_varset{i}.X_h{1};
                end
                tmp=X_h(:);
                X_h_unique=unique(X_h(:));

                cut=sqrt(size(obj.stem_crossval_result.y_back{1},1)*obj.stem_data.T);
                X_h_num=[];
                if length(X_h_unique)/length(tmp)>0.3 
                    flag_mean_h=1;
                    bins=linspace(min(tmp),max(tmp),21);

                    for i=1:(length(bins)-1)
                        bin_h=tmp(tmp>=bins(i)&tmp<bins(i+1));
                        if length(bin_h)>cut
                            X_h_num=cat(1,X_h_num,length(bin_h));
                        end
                    end
                else
                    flag_mean_h=0;
                    for i=1:length(X_h_unique)
                        if sum(tmp==X_h_unique(i))>cut
                            X_h_num=cat(1,X_h_num,sum(tmp==X_h_unique(i)));
                        end
                    end
                end

                % plot of mse_h and mse_t    
                figure 
                if vertical
                    subplot(1,2,1)
                    if length(xvalresult.cv_h)>50
                        plot(xvalresult.cv_mse_h,xvalresult.cv_h,'Color','r' );
                    else
                        if max(X_h_num(:))==min(X_h_num(:))
                             plot(xvalresult.cv_mse_h,xvalresult.cv_h,...
                                'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                                'MarkerSize',6);
                        else
                            scatter(xvalresult.cv_mse_h,xvalresult.cv_h,50,X_h_num,'filled',...
                                'MarkerEdgeColor','k');
                            cl = colorbar;
                            cl.Limits=[min(X_h_num(:)) max(X_h_num(:))];
                        end
                    end
                    if flag_mean_h
                        xlabel('CV-MSE($\bar{h}$)','FontSize',14,'interpreter','latex');
                    else
                        xlabel('CV-MSE(h)','FontSize',14,'interpreter','latex');
                    end
                    set(gca,'YDir','reverse','FontSize',14)
                    if not(isempty(obj.stem_data.stem_varset_p.X_h_unit))
                        ylabel([obj.stem_data.stem_varset_p.X_h_name,' [$',obj.stem_data.stem_varset_p.X_h_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                    else
                        ylabel([obj.stem_data.stem_varset_p.X_h_name],...
                        'FontSize',14,'interpreter','latex');
                    end

                    ylim([xvalresult.cv_h(1),xvalresult.cv_h(end)])
                    grid on
                    box on
                else
                    subplot(1,2,1);
                    if length(xvalresult.cv_h)>50
                        plot(xvalresult.cv_h,xvalresult.cv_mse_h,'Color','r' );
                    else
                        if max(X_h_num(:))==min(X_h_num(:))
                            plot(xvalresult.cv_h,xvalresult.cv_mse_h,...
                                'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                                'MarkerSize',6);
                        else
                            scatter(xvalresult.cv_h,xvalresult.cv_mse_h,50,X_h_num,'filled',...
                                'MarkerEdgeColor','k');
                            cl = colorbar;
                            cl.Limits=[min(X_h_num(:)) max(X_h_num(:))];
                        end
                    end
                    if flag_mean_h
                        ylabel('CV-MSE($\bar{h}$)','FontSize',14,'interpreter','latex');
                    else
                        ylabel('CV-MSE(h)','FontSize',14,'interpreter','latex');
                    end
                    set(gca,'FontSize',14)
                    if not(isempty(obj.stem_data.stem_varset_p.X_h_unit))
                        xlabel([obj.stem_data.stem_varset_p.X_h_name,' [$',obj.stem_data.stem_varset_p.X_h_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                    else
                        xlabel([obj.stem_data.stem_varset_p.X_h_name],...
                        'FontSize',14,'interpreter','latex');
                    end
                    xlim([xvalresult.cv_h(1),xvalresult.cv_h(end)])
                    grid on
                    box on
                end


                subplot(1,2,2);
                if length(xvalresult.cv_mse_t)>50
                    plot(1:length(xvalresult.cv_mse_t), xvalresult.cv_mse_t,...
                        'Color','r' )
                else
                    plot(1:length(xvalresult.cv_mse_t), xvalresult.cv_mse_t,...
                        'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                        'MarkerSize',6)
                end
                ylabel('CV-MSE(t)','FontSize',14,'interpreter','latex');
                set(gca,'FontSize',14)
                if not(isempty(obj.stem_data.stem_varset_p.T_unit))
                    xlabel(['t',' [$',obj.stem_data.stem_varset_p.T_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                else
                    xlabel('t','FontSize',14,'interpreter','latex');
                end

                xlim([1,length(xvalresult.cv_mse_t)])
                grid on
                box on

                % plot of r2_h and r2_t
                figure 
                if vertical
                    subplot(1,2,1)
                    if length(xvalresult.cv_h)>50
                        plot(xvalresult.cv_R2_h,xvalresult.cv_h,'Color','r' );
                    else
                        if max(X_h_num(:))==min(X_h_num(:))
                            plot(xvalresult.cv_R2_h,xvalresult.cv_h,...
                                'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                                'MarkerSize',6);
                        else
                            scatter(xvalresult.cv_R2_h,xvalresult.cv_h,50,X_h_num,'filled',...
                                'MarkerEdgeColor','k');
                            cl = colorbar;
                            cl.Limits=[min(X_h_num(:)) max(X_h_num(:))];
                        end
                    end
                    if flag_mean_h
                        xlabel('CV-R2($\bar{h}$)','FontSize',14,'interpreter','latex');
                    else
                        xlabel('CV-R2(h)','FontSize',14,'interpreter','latex');
                    end
                    set(gca,'YDir','reverse','FontSize',14)
                    if not(isempty(obj.stem_data.stem_varset_p.X_h_unit))
                        ylabel([obj.stem_data.stem_varset_p.X_h_name,' [$',obj.stem_data.stem_varset_p.X_h_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                    else
                        ylabel([obj.stem_data.stem_varset_p.X_h_name],...
                        'FontSize',14,'interpreter','latex');
                    end
                    ylim([xvalresult.cv_h(1),xvalresult.cv_h(end)])
                    grid on
                    box on
                else 
                    subplot(1,2,1)
                    if length(xvalresult.cv_h)>50
                        plot(xvalresult.cv_h,xvalresult.cv_R2_h,'Color','r' );
                    else
                        if max(X_h_num(:))==min(X_h_num(:))
                             plot(xvalresult.cv_h,xvalresult.cv_R2_h,...
                            'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                            'MarkerSize',6);
                        else
                            scatter(xvalresult.cv_h,xvalresult.cv_R2_h,50,X_h_num,'filled',...
                                'MarkerEdgeColor','k');
                            cl = colorbar;
                            cl.Limits=[min(X_h_num(:)) max(X_h_num(:))];
                        end
                    end
                    if flag_mean_h
                        ylabel('CV-R2($\bar{h}$)','FontSize',14,'interpreter','latex');
                    else
                        ylabel('CV-R2(h)','FontSize',14,'interpreter','latex');
                    end
                    set(gca,'FontSize',14)
                    if not(isempty(obj.stem_data.stem_varset_p.X_h_unit))
                        xlabel([obj.stem_data.stem_varset_p.X_h_name,' [$',obj.stem_data.stem_varset_p.X_h_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                    else
                        xlabel([obj.stem_data.stem_varset_p.X_h_name],...
                        'FontSize',14,'interpreter','latex');
                    end
                    xlim([xvalresult.cv_h(1),xvalresult.cv_h(end)])
                    grid on
                    box on
                end

                subplot(1,2,2)
                if length(xvalresult.cv_R2_t)>50
                    plot(1:length(xvalresult.cv_R2_t), xvalresult.cv_R2_t,...
                        'Color','r')
                else
                    plot(1:length(xvalresult.cv_R2_t), xvalresult.cv_R2_t,...
                        'Color','r' ,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                        'MarkerSize',6)
                end
                ylabel('CV-R2(t)','FontSize',14,'interpreter','latex');
                set(gca,'FontSize',14)
                if not(isempty(obj.stem_data.stem_varset_p.T_unit))
                    xlabel(['t',' [$',obj.stem_data.stem_varset_p.T_unit,'$]'],...
                        'FontSize',14,'interpreter','latex');
                else
                    xlabel('t','FontSize',14,'interpreter','latex');
                end
                xlim([1,length(xvalresult.cv_mse_t)])
                grid on
                box on
            end 
            
        end
        
        function plot_profile(obj,lat0,lon0,t_start,t_end)
            %DESCRIPTION: the cross-validation results when modeltype is f-HDGM
            %
            %INPUT
            %obj      - [stem_model object] (1x1) the stem_model object
            %lat0     - [double]            (1x1) the latitude of the spatial site
            %lon0     - [double]            (1x1) the lontitude of the spatial site
            %<t_start>- [double]            (1x1) the start of time step
            %<t_end>  - [double]            (1x1) the end of time step
            %
            %OUTPUT
            %
            %None: the cross-validation results are shown in the figure window
            
            if not(obj.stem_par.stem_modeltype.is('f-HDGM'))
                error('The plot_xval method is only used when model type is f-HDGM')
            end
            
            if nargin<3
                error('not enougth arguments')
            end
            
            if nargin<4
                t_start=[];
            end
            
            if nargin<5
                t_end=[];
            end
            
            [~,idx]=min(abs(obj.stem_data.stem_gridlist_p.grid{1, 1}.coordinate(:,1)-lat0)+...
                abs(obj.stem_data.stem_gridlist_p.grid{1, 1}.coordinate(:,2)-lon0));
            
            lat=obj.stem_data.stem_gridlist_p.grid{1, 1}.coordinate(idx,1);
            lon=obj.stem_data.stem_gridlist_p.grid{1, 1}.coordinate(idx,2);
            
            y_values=[];
            for i=1:length(obj.stem_data.stem_varset_p.Y)
                y_values=cat(1,y_values,obj.stem_data.stem_varset_p.Y{i}(idx,:));
            end
            
            h_values=[];
            for i=1:length(obj.stem_data.stem_varset_p.X_h)
                h_values=cat(1,h_values,obj.stem_data.stem_varset_p.X_h{i}(idx,:));
            end
            
            if isempty(t_start)
                t_values=repmat(1:size(h_values,2),size(h_values,1),1);
                h_values=h_values(:,1:size(h_values,2));
                y_values=y_values(:,1:size(h_values,2));
            elseif isempty(t_end)
                t_values=repmat(t_start:size(h_values,2),size(h_values,1),1);
                h_values=h_values(:,t_start:size(h_values,2));
                y_values=y_values(:,t_start:size(h_values,2));
            else
                t_values=repmat(t_start:t_end,size(h_values,1),1);
                h_values=h_values(:,t_start:t_end);
                y_values=y_values(:,t_start:t_end);
            end
            
            if isempty(obj.stem_data.stem_varset_p.Y_unit)
                disp('Warning: the unit of Y variable is missing.')
                flag_unit=0;
            else
                flag_unit=1;
            end
            figure
            subplot(1,6,1:3);
            scatter(t_values(:),h_values(:),60,y_values(:),'filled',...
                'MarkerEdgeColor','k');
            cl = colorbar;
            cl.Limits=[min(y_values(:)) max(y_values(:))];
            if flag_unit
                cl.Label.String = [obj.stem_data.stem_varset_p.Y_name{1},' [$',obj.stem_data.stem_varset_p.Y_unit,'$]']; 
                cl.Label.Interpreter='latex';
                set(cl,'FontSize',14);
                if isempty(obj.stem_data.stem_varset_p.T_unit)
                    xlabel('t','FontSize',14,'Interpreter','latex');
                else
                    xlabel(['t',' [$',obj.stem_data.stem_varset_p.T_unit,'$]'],'FontSize',14,'Interpreter','latex');
                end
                if isempty(obj.stem_data.stem_varset_p.X_h_unit)
                    ylabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14,'Interpreter','latex'); 
                else
                    ylabel([obj.stem_data.stem_varset_p.X_h_name,' [$',obj.stem_data.stem_varset_p.X_h_unit,'$]'],'FontSize',14,'Interpreter','latex'); 
                end
            else
                cl.Label.String = [obj.stem_data.stem_varset_p.Y_name{1}];
                set(cl,'FontSize',14);
                xlabel('t','FontSize',14);
                ylabel(obj.stem_data.stem_varset_p.X_h_name,'FontSize',14);
            end
            set(gca,'Ydir','reverse','FontSize',14);
            ylim([min(h_values(:)), max(h_values(:))]);
            xlim([min(t_values(:))-1, max(t_values(:))+1]);
            title({'Nearest location to coordinates provided is',['Lat ', num2str(lat),' and Lon ', num2str(lon)]},'FontSize',14);
            box on
            grid on
            
            subplot(1,6,4.5:5.8);
            plot(y_values,h_values,'.-')
            if flag_unit
                xlabel([obj.stem_data.stem_varset_p.Y_name{1},' [$',...
                    obj.stem_data.stem_varset_p.Y_unit,'$]'],'FontSize',14,'Interpreter','latex');
            else
                xlabel(obj.stem_data.stem_varset_p.Y_name{1},'FontSize',14);
            end
            set(gca,'Ydir','reverse','FontSize',14);
            ylim([min(h_values(:)), max(h_values(:))]);
            xlim([min(y_values(:))-1, max(y_values(:))+1]);
            box on
            grid on
        end
        
        function beta_Chi2_test(obj)
            %DESCRIPTION: print the Chi2 test for significance of covariates when modeltype is f-HDGM
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: the test result is printed in the command window
            
            if not(obj.stem_par.stem_modeltype.is('f-HDGM')&&obj.stem_par.flag_beta_spline==1)
                error('The beta_Chi2_test function is only used when model type is f-HDGM and beta is spline coefficient')
            end
            if isempty(obj.stem_par.varcov)
                error('We can not get the Chi2 statistic if variance covariance matrix is not estimated')
            end
            
            k = obj.stem_par.k_beta;
            output=cell(1+length(obj.stem_data.stem_varset_p.X_beta_name{1,1}),3);
            output{1,1}='Covariate';
            output{1,2}='Chi2';
            output{1,3}='P-Value';
            for i = 1:length(obj.stem_data.stem_varset_p.X_beta_name{1,1})
                output{i+1,1} = obj.stem_data.stem_varset_p.X_beta_name{1,1}{i};            
                c_hat = obj.stem_par.beta((i-1)*k+(1:k));
                sigma_c_hat=obj.stem_par.varcov((i-1)*k+(1:k),(i-1)*k+(1:k));
                chi2_hat=c_hat'*(inv(sigma_c_hat))*c_hat;
                output{i+1,2}=num2str(chi2_hat);
                output{i+1,3}=num2str(1-chi2cdf(chi2_hat,k));
            end
            disp(output);   
        end  
        
        function print2xlsx(obj,filename)
            %DESCRIPTION: print the information on the estimation results
            %to a xlsx file
            %
            %INPUT
            %obj        - [stem_model object]   (1x1) the stem_data object
            %filename   - [string]              (1x1) the name of the xlsx file
            %
            %OUTPUT
            %
            %none: the information is printed in the xlsx file
            sheet='D-STEM';
            r=1;
            c=1;
            if obj.estimated
                if obj.tapering
                    string='Tapering is enabled'; 
                    xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                    r=r+1;
                    if not(isempty(obj.tapering_p))
                        string=['Point data tapering: ',num2str(obj.stem_data.stem_gridlist_p.tap),' km'];
                    else
                        string='Tapering is NOT enabled on point data';
                    end
                    xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                    r=r+1;
                    if not(isempty(obj.tapering_b))
                        string=['Pixel data tapering: ',num2str(obj.stem_data.stem_gridlist_b.tap),' km'];
                    else
                        string='Tapering is NOT enabled on pixel data';
                    end
                    xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                    r=r+1;
                else
                    string='Tapering is not enabled';
                    xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                    r=r+1;
                end
                r=r+1;
                if not(isempty(obj.stem_EM_result.logL))
                    string=['Observed data log-likelihood: ',num2str(obj.stem_EM_result.logL)];
                else
                    string='Observed data log-likelihood: not computed';
                end
                xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                r=r+2;
                string='Beta of point variables';
                xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                r=r+1;
                counter=1;
                counter_varcov=1;
                if not(isempty(obj.stem_par.beta))
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output=cell(length(obj.stem_data.stem_varset_p.X_beta_name{i})+2,4);
                        output{1,1}=['Beta of ',obj.stem_data.stem_varset_p.Y_name{i}];
                        output{2,1}='Loading';
                        output{2,2}='Value';
                        output{2,3}='Std';
                        output{2,4}='|t|';
                        for j=1:length(obj.stem_data.stem_varset_p.X_beta_name{i})
                            output{j+2,1}=obj.stem_data.stem_varset_p.X_beta_name{i}{j};
                            output{j+2,2}=num2str(obj.stem_par.beta(counter));
                            if not(isempty(obj.stem_par.varcov))
                                output{j+2,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                            else
                                output{j+2,3}='Not computed';
                            end
                            if not(isempty(obj.stem_par.varcov))
                                output{j+2,4}=num2str(abs(obj.stem_par.beta(counter)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))));
                            else
                                output{j+2,4}='Not computed';
                            end
                            counter=counter+1;
                            counter_varcov=counter_varcov+1;
                        end
                        xlswrite(filename,output,sheet,stem_misc.rc2xls(r,(i-1)*4+1));
                    end
                    r=r+size(output,1)+2;
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_beta))
                            string='Beta of pixel variables';
                            xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                            r=r+1;
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                output=cell(length(obj.stem_data.stem_varset_b.X_beta_name{i})+1,4);
                                output{1,1}=['Beta of ',obj.stem_data.stem_varset_b.Y_name{i}];
                                output{2,1}='Loading coefficient';
                                output{2,2}='Value';
                                output{2,3}='Std';
                                output{2,4}='|t|';
                                for j=1:length(obj.stem_data.stem_varset_b.X_beta_name{i})
                                    output{j+2,1}=obj.stem_data.stem_varset_b.X_beta_name{i}{j};
                                    output{j+2,2}=num2str(obj.stem_par.beta(counter));
                                    if not(isempty(obj.stem_par.varcov))
                                        output{j+2,3}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                                    else
                                        output{j+2,3}='Not computed';
                                    end
                                    if not(isempty(obj.stem_par.varcov))
                                        output{j+2,4}=num2str(abs(obj.stem_par.beta(counter)/sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov))));
                                    else
                                        output{j+2,4}='Not computed';
                                    end
                                    counter=counter+1;
                                    counter_varcov=counter_varcov+1;
                                end
                                xlswrite(filename,output,sheet,stem_misc.rc2xls(r,(i-1)*4+1));
                            end
                            r=r+size(output,1)+2;
                        end
                    end
                end

                output=cell(4,obj.stem_data.stem_varset_p.nvar+1);
                output{1,1}='sigma_eps point variables';
                output{2,1}='Variable';
                output{3,1}='Value';
                output{4,1}='Std';
                for i=1:obj.stem_data.stem_varset_p.nvar
                    output{2,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                    output{3,i+1}=num2str(obj.stem_par.sigma_eps(i,i));
                    if not(isempty(obj.stem_par.varcov))
                        output{4,i+1}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                    else
                        output{4,i+1}='Not computed';
                    end
                    counter=counter+1;
                    counter_varcov=counter_varcov+1;
                end
                xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c));
                r=r+size(output,1)+2;

                if not(isempty(obj.stem_data.stem_varset_b))
                    output=cell(4,obj.stem_data.stem_varset_b.nvar+1);
                    output{1,1}='sigma_eps pixel variables';
                    output{2,1}='Variable';
                    output{3,1}='Value';
                    output{4,1}='Std';
                    delta=obj.stem_data.stem_varset_p.nvar;
                    for i=1:obj.stem_data.stem_varset_b.nvar
                        output{2,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                        output{3,i+1}=num2str(obj.stem_par.sigma_eps(i+delta,i+delta));
                        if not(isempty(obj.stem_par.varcov))
                            output{4,i+1}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                        else
                            output{4,i+1}='Not computed';
                        end
                        counter=counter+1;
                        counter_varcov=counter_varcov+1;
                    end
                    xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c));
                    r=r+size(output,1)+2;
                end
                
                if not(isempty(obj.stem_data.stem_varset_b))
                   
                    
                    output=[];
                    if obj.stem_par.stem_par_constraints.pixel_correlated
                        string='Pixel variables are cross-correlated';
                        xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                        r=r+1;
                        output{1,1}='theta_b';
                        output{2,1}='Value [km]';
                        output{3,1}='Std [km]';
                        output{2,2}=num2str(obj.stem_par.theta_b(1));
                        if not(isempty(obj.stem_par.varcov))
                            output{3,2}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                        else
                            output{3,2}='Not computed';
                        end
                        counter=counter+1;
                        counter_varcov=counter_varcov+1;
                        xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c));
                        r=r+size(output,1)+2;
                        
                        string='v_b';
                        xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                        r=r+1;
                        
                        output=cell(obj.stem_data.stem_varset_b.nvar+1,obj.stem_data.stem_varset_b.nvar+1);
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,i+1}=num2str(1);
                        end
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            for j=i+1:obj.stem_data.stem_varset_b.nvar
                                output{i+1,j+1}=num2str(obj.stem_par.v_b(i,j));
                                counter=counter+1;
                            end
                        end
                        xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c));
                        if not(isempty(obj.stem_par.varcov))
                            string='v_b std';
                            xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c+obj.stem_data.stem_varset_b.nvar+2));
                            output=cell(obj.stem_data.stem_varset_b.nvar+1,obj.stem_data.stem_varset_b.nvar+1);
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                output{1,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                                output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                                output{i+1,i+1}=num2str(1);
                            end
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                for j=i+1:obj.stem_data.stem_varset_b.nvar
                                    output{i+1,j+1}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                                    counter_varcov=counter_varcov+1;
                                end
                            end
                            xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c+obj.stem_data.stem_varset_b.nvar+2));
                        end
                        r=r+size(output,1)+2;
                    else
                        string='Pixel variables are not cross-correlated';
                        xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
                        r=r+1;

                        output=cell(4,obj.stem_data.stem_varset_b.nvar+1);
                        output{1,1}='Theta_b';
                        output{2,1}='Variable';
                        output{3,1}='Value [km]';
                        output{4,1}='Std [km]';
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{2,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{3,i+1}=num2str(obj.stem_par.theta_b(i));
                            if not(isempty(obj.stem_par.varcov))
                                output{4,i+1}=num2str(sqrt(obj.stem_par.varcov(counter_varcov,counter_varcov)));
                            else
                                output{4,i+1}='Not computed';
                            end
                            counter=counter+1;
                            counter_varcov=counter_varcov+1;
                        end
                        xlswrite(filename,output,sheet,stem_misc.rc2xls(r,c));
                        r=r+size(output,1)+2;
                    end
                end
                
                if obj.stem_par.k>0
                    string=['The model has ',num2str(obj.stem_par.k),' fine-scale coregionalization component(s) w_p'];
                    xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));                    
                    r=r+1;
                 end


            else
                string='The model is not estimated';
                xlswrite(filename,{string},sheet,stem_misc.rc2xls(r,c));
            end
        end

        function all_par = par_vec(obj)
            %DESCRIPTION: vectorize the model parameters (only the estimated parameters with the exclusion of structural zeroes and repeated elements of symmetric matrices)
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %all_par - [double]            (Hx1) the vectorized the model parameters
                  
            all_par=[];
            all_par=[all_par; obj.stem_par.beta];
            if size(obj.stem_par.sigma_eps,3)==1
                if obj.stem_par.flag_sigma_eps_spline==1 
                    all_par=[all_par; obj.stem_par.sigma_eps];
                else
                    all_par=[all_par; diag(obj.stem_par.sigma_eps)];
                end 
            end

           
            all_par=[all_par; obj.stem_par.theta_b(:)];

            if obj.stem_par.stem_par_constraints.pixel_correlated
                all_par=[all_par; stem_misc.from_upper_triangular_to_vector(obj.stem_par.v_b)];
            end

            
            all_par=[all_par;obj.stem_par.theta_p(:)];
            for i=1:size(obj.stem_par.v_p,3)
                all_par=cat(1,all_par,stem_misc.from_upper_triangular_to_vector(obj.stem_par.v_p(:,:,i)));
            end
            
            if obj.stem_par.p>0
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    if not(obj.stem_par.stem_par_constraints.time_diagonal)
                        all_par=[all_par; obj.stem_par.G(:)];
                        all_par=[all_par; stem_misc.triuv(obj.stem_par.sigma_eta)];
                    else
                        all_par=[all_par; diag(obj.stem_par.G)];
                        all_par=[all_par; diag(obj.stem_par.sigma_eta)];
                    end
                else
                    all_par=[all_par; diag(obj.stem_par.G)];
                end
            end
            
            if obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})      
                all_par=[all_par;obj.stem_par.theta_z(:)];
                all_par=[all_par;diag(obj.stem_par.v_z)];
                %all_par=[all_par;stem_misc.from_upper_triangular_to_vector(obj.stem_par.v_z)];
            end
        end
        
        function print_par(obj)
            %DESCRIPTION: print the stem_par object
            %
            %INPUT
            %obj     - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %none: the stem_par object is printed in the command window 
            
            disp('****************');
            disp('PARAMETER VALUES');
            disp('****************');
            
            
            if not(isempty(obj.stem_par.beta))
                disp(['beta: ',num2str(obj.stem_par.beta')]);
            end
            if obj.stem_par.flag_sigma_eps_spline==1 
                disp(['sigma_eps: ',num2str(obj.stem_par.sigma_eps')]);
            else
                if size(obj.stem_par.sigma_eps,3)==1
                    disp(['sigma_eps: ',num2str(diag(obj.stem_par.sigma_eps)')]);
                else
                    disp(['Time avg. sigma_eps: ',num2str(diag(mean(obj.stem_par.sigma_eps,3))')]);
                end
            end
            
            if not(isempty(obj.stem_par.alpha_bp))
                disp(['alpha_bp: ',num2str(obj.stem_par.alpha_bp')]);
            end
            if not(isempty(obj.stem_par.theta_b))
                disp(['theta_b: ',num2str(obj.stem_par.theta_b(:)')]);
            end
            if not(isempty(obj.stem_par.theta_p))
                disp(['theta_p: ',num2str(obj.stem_par.theta_p(:)')]);
            end
            if not(isempty(obj.stem_par.v_b))
                if (obj.stem_par.stem_par_constraints.pixel_correlated)&&(obj.stem_par.q>1)
                    disp('v_b: ');
                    disp(obj.stem_par.v_b);
                end
            end
            if not(isempty(obj.stem_par.v_p))
                for i=1:min(size(obj.stem_par.v_p,3),obj.stem_par.k)
                    disp(['v_p',num2str(i),': ']);
                    disp(obj.stem_par.v_p(:,:,i));
                end
            end
            
            if obj.stem_par.p>0
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    if obj.stem_par.stem_par_constraints.time_diagonal
                        disp(['diag_G: ',num2str(diag(obj.stem_par.G)')]);
                        disp(['sigma_eta: ',num2str(diag(obj.stem_par.sigma_eta)')]);
                    else
                        disp('G: ');
                        disp(obj.stem_par.G);
                        disp('sigma_eta: ');
                        disp(obj.stem_par.sigma_eta);
                    end
                    if not(isempty(obj.stem_par.lambda))
                        disp(['lambda: ',num2str(obj.stem_par.lambda)]);
                    end
                else
                    disp(['diag_G: ',num2str(diag(obj.stem_par.G)')]);
                    disp(['theta_z: ',num2str(obj.stem_par.theta_z(:)')]);
                    disp('v_z: ');
                    disp(obj.stem_par.v_z);
                end
            end
        end     
        
        function set_system_size(obj)
            %DESCRIPTION: evaluate the minimum N after which it is faster to evaluate only the diagonal of a matrix product instead of the full matrix
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: the system_size property is updated
            
            dim=obj.N;
            if dim>200
                n_blocks=[1 2 3 5 10 15 25 50 100 150];
                n_blocks=[1];
                times=zeros(length(n_blocks),1);
                mat1=randn(dim);
                mat2=randn(dim);
                res=zeros(dim,1);
                last_step=0;
                counter=1;
                for i=1:length(n_blocks)
                    step=ceil(dim/n_blocks(i));
                    if not(step==last_step)
                        step_vec(i)=step;
                        
                        blocks=0:step:dim;
                        if not(blocks(end)==dim)
                            blocks=[blocks dim];
                        end
                        tic;
                        for z=1:max(1,floor(10000/dim))
                            for j=1:length(blocks)-1
                                res(blocks(j)+1:blocks(j+1),1)=diag(mat1(blocks(j)+1:blocks(j+1),:)*mat2(blocks(j)+1:blocks(j+1),:)');
                            end
                        end
                        t=toc;
                        times(i)=t;
                        last_step=step;
                        counter=counter+1;
                    end
                end
                [~,idx]=min(times(1:counter-1));
                if idx>1
                    obj.product_step=step_vec(idx);
                else
                    obj.product_step=-1;
                end
            else
                obj.product_step=-1;
            end
        end
        
        function set_logL(obj)
            %DESCRIPTION: evaluate the observed data log-likelihood
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: the logL property of the object stem_EM_result is updated
            disp('Log-Likelihood computation...');
            st_kalman=stem_kalman(obj);
            
            [st_kalmansmoother_result,~,~,~,~,~,~] = st_kalman.smoother(1);
            obj.stem_EM_result.logL=st_kalmansmoother_result.logL;
            disp('Log-Likelihood computation ended.');
        end
        
        function set_varcov(obj,exit_toll)
            %DESCRIPTION: evaluate the variance-covariance matrix of the model parameters
            %
            %INPUT
            %
            %obj       - [stem_model object] (1x1) the stem_model object
            %exit_toll - [double>0]          (1x1) exit tolerance for approximated variance-covariance matrix
            %
            %OUTPUT
            %
            %none: the varcov property of the object stem_EM_result is updated
            
            %parameter order: beta,sigma_eps,alpha_bp,theta_b,v_b,alpha_p,theta_p,v_p,G,sigma_eta,[alpha_z,theta_z,v_z]
            
            if obj.estimated==0
                error('The model has not been estimated yet');
            end
            
            if nargin>1
                if exit_toll<=0
                    error('The input exit_toll must be > 0');
                end
            else
                exit_toll=0;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  dimension estimation   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            N=obj.N;
            T=obj.T;
            dim=obj.dim;
            
            data=obj.stem_data;
            par=obj.stem_par;
            q=par.q;
            p=par.p;
            k=par.k;
            
            if obj.stem_data.stem_modeltype.is('f-HDGM')
                Np=p*dim(1);
            else
                Np=obj.Np;
            end
            
            if obj.stem_data.stem_modeltype.is('Emulator')
                d=length(obj.DistMat_p);
            else
                d=0;
            end
            
            sigma_eps_estimated=1-par.stem_par_constraints.sigma_eps_fixed;

            if p>0
                temp=abs(par.G-eye(size(par.G)));
                if sum(temp(:))==0
                    random_walk=1;
                else
                    random_walk=0;
                end
            end
            if k>0
                if k>1&&size(par.v_p,3)==1
                    shared_v=1;
                else
                    shared_v=0;
                end
            end
            
            n_beta=0;
            n_bp_theta=0;
            n_bp_v=0;
            n_bp_alpha=0;
            n_p_theta=0;
            n_p_v=0;
            n_time_s2e=0;

            if not(isempty(data.X_bp))
                n_bp_alpha=q;
                if par.stem_par_constraints.pixel_correlated
                    n_bp_v=q*(q+1)/2;
                    if not(strcmp(par.correlation_type,'expsphere'))
                        n_bp_theta=1;
                    else
                        n_bp_theta=2;
                    end
                else
                    if not(strcmp(par.correlation_type,'expsphere'))
                        n_bp_theta=q;
                    else
                        n_bp_theta=2*q;
                    end
                    n_bp_v=q;
                end
            end
            
            if not(isempty(data.X_beta))
                n_beta=length(par.beta);
            end
            
            if p>0
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    if par.stem_par_constraints.time_diagonal
                        if random_walk
                            n_time_G=0;
                        else
                            n_time_G=p;
                        end
                        n_time_s2e=p;
                    else
                        n_time_G=p^2;
                        n_time_s2e=(p+1)*p/2;
                    end
                    n_time=n_time_G+n_time_s2e;
                    
                else
                    if random_walk
                        n_time_G=0;
                    else
                        n_time_G=p;
                    end
                   
                    if not(strcmp(par.correlation_type,'expsphere'))
                        if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                            n_z_theta=1;
                        else
                            n_z_theta=p;
                        end
                    else
                        if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                            n_z_theta=2;
                        else
                            n_z_theta=p*2;
                        end
                    end
                    if obj.stem_data.stem_modeltype.is('f-HDGM')
                        n_z_v=p;
                    else
                        n_z_v=(p+1)*p/2;
                    end
                    n_time_s2e=0;
                    n_time=n_time_G+n_z_theta+n_z_v;
                end
            else
                n_time=0;
            end
            
            if not(isempty(data.X_p))
                if shared_v
                    n_p_v=q*(q+1)/2;
                else
                    n_p_v=q*(q+1)/2*k;
                end
                if not(obj.stem_data.stem_modeltype.is('Emulator'))
                    if not(strcmp(par.correlation_type,'expsphere'))
                        n_p_theta=k;
                    else
                        n_p_theta=2*k;
                    end
                else
                    n_p_theta=d;
                end
            end
            if sigma_eps_estimated
                n_eps=size(par.sigma_eps,1)*size(par.sigma_eps,3);
            else
                n_eps=0;
            end
            
            n_psi=n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_theta+n_p_v+n_time;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  extraction of the useful data   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if p>0
                %kalman filter
                st_kalman=stem_kalman(obj);
                compute_logL=0;
                enable_varcov_computation=1;
                [st_kalmanfilter_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_eta,G_tilde_diag,sigma_geo,aj_bp,M] = st_kalman.filter(compute_logL,enable_varcov_computation);
                rr=size(sigma_Z,1);
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    G=par.G;
                else
                    G=sparse(1:length(G_tilde_diag),1:length(G_tilde_diag),G_tilde_diag,length(G_tilde_diag),length(G_tilde_diag));
                end
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,sigma_eta,~,aj_bp,M] = obj.get_sigma();
                st_kalmanfilter_result=stem_kalmanfilter_result([],[],[],[],[],[],[]);
                rr=0;
            end
            J=st_kalmanfilter_result.J;
            zk_f=st_kalmanfilter_result.zk_f;
            Pk_f=st_kalmanfilter_result.Pk_f;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  derivatives allocation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Derivatives allocation...');
            d_e_Lt=cell(n_psi,1);
            if not(isempty(data.X_z))
                d_Z=cell(n_psi,1);
                d_J_Lt=cell(n_psi,1);
                d_P=cell(n_psi,1);
            end
            d_St_Lt=cell(n_psi,1);
            
            
            d_Sgeo=cell(n_psi,size(par.sigma_eps,3));
            
            if not(isempty(data.X_beta))
                d_beta=cell(n_psi,1);
            end
            if sigma_eps_estimated
                d_Seps=cell(n_psi,size(par.sigma_eps,3));
            end
            if not(isempty(data.X_bp))
                d_alpha_bp=cell(n_psi,n_bp_alpha);
                d_M_Sb=cell(n_psi,1);
            end
            if not(isempty(data.X_p))
                d_Sp=cell(n_psi,k);
            end
            if not(isempty(data.X_z))
                d_G=cell(n_psi,1);
                d_s2e=cell(n_psi,1);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:n_psi
                if not(isempty(data.X_beta))
                    d_beta{i}=sparse(n_beta,1);
                end
                if sigma_eps_estimated 
                    if par.flag_sigma_eps_spline==1
                        d_Seps{i} = sparse(n_eps,1);
                    else
                        for t=1:size(par.sigma_eps,3)
                            d_Seps{i,t} = sparse(N,N);
                        end
                    end
                end
                if not(isempty(data.X_bp))
                    for j=1:n_bp_alpha
                        d_alpha_bp{i,j}=sparse(N,1);
                    end
                    d_M_Sb{i}=sparse(N,N);
                end
                if not(isempty(data.X_p))
                    for j=1:k
                        d_Sp{i,j}=sparse(Np,Np);
                    end
                end
                if not(isempty(data.X_z))
                    if obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                        if not(isempty(obj.DistMat_p))
                            if not(iscell(obj.DistMat_p))
                                NN=size(obj.DistMat_p,1);
                            else
                                NN=size(obj.DistMat_p{1},1);
                            end
                        end
                        if not(isempty(obj.DistMat_z))
                            if not(iscell(obj.DistMat_z))
                                NN=size(obj.DistMat_z,1);
                            else
                                NN=size(obj.DistMat_z{1},1);
                            end
                        end
                        
                    end
                    d_G{i}=sparse(rr,rr);
                    d_s2e{i}=sparse(rr,rr);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  initialization of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Derivatives initialization...');
           
            blocks=[0 cumsum(dim)];
            
            %d_beta
            if not(isempty(data.X_beta))
                for i=1:n_beta
                    d_beta{i}(i,1)=1;
                end
            end
            
            %d_Seps
            if sigma_eps_estimated
                if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                    if size(sigma_eps,3)==1
                        for i=1:n_eps
                            Id=blocks(i)+1:blocks(i+1);
                            d_Seps{n_beta+i}=sparse(Id,Id,ones(length(Id),1),N,N);
                        end
                    else
                        for t=1:T
                            for i=1:q
                                Id=blocks(i)+1:blocks(i+1);
                                d_Seps{n_beta+(t-1)*q+i,t}=sparse(Id,Id,ones(length(Id),1),N,N);
                            end
                        end
                    end
                else 
                    if par.flag_sigma_eps_spline==1
                        for i=1:n_eps
                            d_Seps{n_beta+i}(i,1) = 1;
                        end
                     else
                        d_Seps{n_beta+1}=speye(N);
                     end
                end
            end
            
            if not(isempty(data.X_bp))
                %d_alpha_bp
                for i=1:n_bp_alpha
                    [~,j_b] = obj.get_jbp(i);
                    d_alpha_bp{n_beta+n_eps+i,i}=j_b;
                end
                
                %d_M_Sb with respect to theta_bp
                if not(strcmp(par.correlation_type,'expsphere'))
                    d=stem_misc.M_apply(obj.DistMat_b,M,'b');
                    result=stem_misc.M_apply(sigma_W_b,M,'b');
                    if not(par.stem_par_constraints.pixel_correlated)&&(q>1)
                        for j=1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            for i=1:2
                                idx=blocks(j+(i-1)*q)+1:blocks(j+(i-1)*q+1);
                                result1=d(idx,idx);
                                result2=result(idx,idx);
                                if strcmp(par.correlation_type,'exponential')
                                    result3=result1/(par.theta_b(j)^2).*result2;
                                elseif strcmp(par.correleation_type,'matern32')
                                    result3=-sqrt(3)*result1/(par.theta_b(j)^2).*result2./(1+sqrt(3).*result1./par.theta_b(j))+result2.*sqrt(3).*result1./(par.theta_b(j)^2);
                                else
                                    result3=(-sqrt(5)*result1/(par.theta_b(j)^2)-10/3*result1.^2/par.theta_b(j)^3).*result2./(1+sqrt(5).*result1./par.theta_b(j)+5/3*result1.^2/par.theta_b(j)^2)+result2.*sqrt(5).*result1./(par.theta_b(j)^2);
                                end
                                
                                L=find(result3);
                                [idx_I,idx_J]=ind2sub(size(result3),L);
                                Id=cat(1,Id,idx_I+blocks(j+(i-1)*2));
                                Jd=cat(1,Jd,idx_J+blocks(j+(i-1)*2));
                                elements=cat(1,elements,result3(L));
                            end
                            
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>90
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+j}=sparse(Id,Jd,elements,N,N);
                            else
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+j}=zeros(N);
                                temp=sub2ind([N,N],Id,Jd);
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+j}(temp)=elements;
                            end
                        end
                    else
                        d_M_Sb{n_beta+n_eps+n_bp_alpha+1}=result.*d/(par.theta_b^2);
                    end
                else
                    for i=1:2
                        d{i}=stem_misc.M_apply(obj.DistMat_b{i},M,'b');
                    end
                    result=stem_misc.M_apply(sigma_W_b,M,'b');
                    if not(par.stem_par_constraints.pixel_correlated)&&(q>1)
                        counter=1;
                        for h=1:2 %loop on the theta parameters of the expsphere correlation function
                            for j=1:q
                                Id=[];
                                Jd=[];
                                elements=[];
                                for i=1:2
                                    idx=blocks(j+(i-1)*q)+1:blocks(j+(i-1)*q+1);
                                    result1=d{h}(idx,idx);
                                    result2=result(idx,idx);
                                    result3=result1/(par.theta_b(h,j)^2).*result2;
                                    
                                    L=find(result3);
                                    [idx_I,idx_J]=ind2sub(size(result3),L);
                                    Id=cat(1,Id,idx_I+blocks(j+(i-1)*2));
                                    Jd=cat(1,Jd,idx_J+blocks(j+(i-1)*2));
                                    elements=cat(1,elements,result3(L));
                                end
                                
                                zero_density=(1-length(Id)/(N^2))*100;
                                if obj.tapering||zero_density>90
                                    d_M_Sb{n_beta+n_eps+n_bp_alpha+counter}=sparse(Id,Jd,elements,N,N);
                                else
                                    d_M_Sb{n_beta+n_eps+n_bp_alpha+counter}=zeros(N);
                                    temp=sub2ind([N,N],Id,Jd);
                                    d_M_Sb{n_beta+n_eps+n_bp_alpha+counter}(temp)=elements;
                                end
                                counter=counter+1;
                            end
                        end
                    else
                        d_M_Sb{n_beta+n_eps+n_bp_alpha+1}=result.*d{1}/(par.theta_b(1)^2);
                        d_M_Sb{n_beta+n_eps+n_bp_alpha+2}=result.*(d{2}.^2)/(par.theta_b(2)^2);
                    end
                end
                
                %d_M_Sb with respect to v_bp
                if par.stem_par_constraints.pixel_correlated
                    z=1;
                    for j=1:q
                        for i=j:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            
                            for h=1:2
                                idx1=blocks(i+(h-1)*q)+1:blocks(i+(h-1)*q+1);
                                idx2=blocks(j+(h-1)*q)+1:blocks(j+(h-1)*q+1);
                                result1=result(idx1,idx2);
                                result2=result1/par.v_b(i,j);
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=cat(1,Id,idx_I+blocks(i+(h-1)*q));
                                Jd=cat(1,Jd,idx_J+blocks(j+(h-1)*q));
                                elements=cat(1,elements,result2(L));
                                if i~=j
                                    Id=cat(1,Id,idx_J+blocks(j+(h-1)*q));
                                    Jd=cat(1,Jd,idx_I+blocks(i+(h-1)*q));
                                    elements=cat(1,elements,result2(L));
                                end
                            end
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>90
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=sparse(Id,Jd,elements,N,N);
                            else
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=zeros(N);
                                temp=sub2ind([N,N],Id,Jd);
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                else
                    z=1;
                    for i=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        
                        for h=1:2
                            idx1=blocks(i+(h-1)*q)+1:blocks(i+(h-1)*q+1);
                            idx2=blocks(i+(h-1)*q)+1:blocks(i+(h-1)*q+1);
                            result1=result(idx1,idx2);
                            result2=result1/par.v_b(i,i);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=cat(1,Id,idx_I+blocks(i+(h-1)*q));
                            Jd=cat(1,Jd,idx_J+blocks(i+(h-1)*q));
                            elements=cat(1,elements,result2(L));
                        end
                        
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>90
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=sparse(Id,Jd,elements,N,N);
                        else
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=zeros(N);
                            temp=sub2ind([N,N],Id,Jd);
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}(temp)=elements;
                        end
                        z=z+1;
                    end
                end
            end
            
            if not(isempty(data.X_p))
                %d_Sp with respect to theta_p
                if not(obj.stem_data.stem_modeltype.is('Emulator'))
                    if not(strcmp(par.correlation_type,'expsphere'))
                        for z=1:k
                            if strcmp(par.correlation_type,'exponential')
                                derivative=sigma_W_p{z}.*obj.DistMat_p/(par.theta_p(z)^2);
                            elseif strcmp(par.correlation_type,'matern32')
                                derivative=-sqrt(3)*obj.DistMat_p/(par.theta_p(z)^2).*exp(-sqrt(3)*obj.DistMat_p/par.theta_p(z))+sigma_W_p{z}.*sqrt(3).*obj.DistMat_p/(par.theta_p(z)^2);
                            else
                                derivative=(-sqrt(5)*obj.DistMat_p/par.theta_p(z)^2-10/3.*obj.DistMat_p.^2/par.theta_p(z)^3).*exp(-sqrt(5)*obj.DistMat_p/par.theta_p(z))+sigma_W_p{z}.*sqrt(5).*obj.DistMat_p/(par.theta_p(z)^2);
                            end
                            if stem_misc.zero_density(derivative)>90
                                derivative=sparse(derivative);
                            end
                            d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z,z}=derivative;
                        end
                    else
                        counter=1;
                        for z=1:k
                            derivative=sigma_W_p{z}.*obj.DistMat_p{1}/(par.theta_p(1,z)^2);
                            if stem_misc.zero_density(derivative)>90
                                derivative=sparse(derivative);
                            end
                            d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+counter,counter}=derivative;
                            counter=counter+1;
                            derivative=sigma_W_p{z}.*obj.DistMat_p{2}/(par.theta_p(2,z)^2);
                            if stem_misc.zero_density(derivative)>90
                                derivative=sparse(derivative);
                            end
                            d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+counter,counter}=derivative;
                            counter=counter+1;
                        end
                    end
                else
                    for z=1:d
                        if strcmp(par.correlation_type,'exponential')
                            derivative=sigma_W_p{1}.*obj.DistMat_p{z}/(par.theta_p(z)^2);
                        elseif strcmp(par.correlation_type,'matern32')
                            derivative=-sqrt(3)*obj.DistMat_p{z}/(par.theta_p(z)^2).*exp(-sqrt(3)*obj.DistMat_p{z}/par.theta_p(z))+sigma_W_p{1}.*sqrt(3).*obj.DistMat_p{z}/(par.theta_p(z)^2);
                        else
                            derivative=(-sqrt(5)*obj.DistMat_p{z}/par.theta_p(z)^2-10/3.*obj.DistMat_p{z}.^2/par.theta_p(z)^3).*exp(-sqrt(5)*obj.DistMat_p{z}/par.theta_p(z))+sigma_W_p{1}.*sqrt(5).*obj.DistMat_p{z}/(par.theta_p(z)^2);
                        end
                        if stem_misc.zero_density(derivative)>90
                            derivative=sparse(derivative);
                        end
                        d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+z}=derivative;
                    end
                end
                
                %d_Sp with respect to v_p
                counter=1;
                for z=1:k
                    for j=1:q
                        for i=j:q    
                            Id=[];
                            Jd=[];
                            elements=[];

                            result1=sigma_W_p{z}(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                            result2=result1/par.v_p(i,j,z);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=cat(1,Id,idx_I+blocks(j));
                            Jd=cat(1,Jd,idx_J+blocks(i));
                            elements=cat(1,elements,result2(L));
                            if i~=j
                                Id=cat(1,Id,idx_J+blocks(i));
                                Jd=cat(1,Jd,idx_I+blocks(j));
                                elements=cat(1,elements,result2(L));
                            end
                            
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>90
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_theta+counter,z}=sparse(Id,Jd,elements,Np,Np);
                            else
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_theta+counter,z}=zeros(Np);
                                temp=sub2ind([Np,Np],Id,Jd);
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_theta+counter,z}(temp)=elements;
                            end
                            counter=counter+1;
                        end
                    end
                end
            end
            
            if not(isempty(data.X_z))
                if not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    %d_G
                    if par.stem_par_constraints.time_diagonal
                        if random_walk==0
                            for i=1:rr
                                d_G{n_psi-n_time+i}(i,i)=1;
                            end
                        end
                    else
                        for i=1:rr^2
                            [a,b]=ind2sub([rr rr],i);
                            d_G{n_psi-n_time+i}(a,b)=1;
                        end
                    end
                    %d_s2e
                    if par.stem_par_constraints.time_diagonal   
                        for i=1:p
                            d_s2e{n_psi-n_time_s2e+i}(i,i)=1;
                        end
                    else
                        j=1;
                        for i=1:(rr*(rr+1))/2
                            temp=zeros(rr);
                            l=1;
                            for h=1:rr
                                for z=h:rr
                                    if i==l
                                        temp(h,z)=1;
                                        temp(z,h)=1;
                                    end
                                    l=l+1;
                                end
                            end
                            d_s2e{n_psi-n_time_s2e+i}=sparse(temp);
                            j=j+1;
                        end
                    end
                else
                    %d_G
                    if random_walk==0
                        dim=obj.stem_data.dim;
                        if p==q
                            blocks=[0 cumsum(dim(1:p))];
                        else
                            blocks=[0 cumsum(ones(1,p)*dim(1))];
                        end
                        for i=1:p
                            j_z=zeros(rr,1);
                            j_z(blocks(i)+1:blocks(i+1))=1;
                            d_G{n_psi-n_time+i}=sparse(1:rr,1:rr,j_z,rr,rr);
                        end
                    end

                    %d_s2e
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        if strcmp(par.correlation_type,'exponential')
                            d_s2e{n_psi-n_time+p+1}=sigma_eta.*obj.DistMat_p/(par.theta_z^2);
                        elseif strcmp(par.correlation_type,'matern32')
                            d_s2e{n_psi-n_time+p+1}=-sqrt(3)*obj.DistMat_p/(par.theta_z^2).*exp(-sqrt(3)*obj.DistMat_p/par.theta_z)+sigma_eta.*sqrt(3).*obj.DistMat_p/(par.theta_z^2);
                        elseif strcmp(par.correlation_type,'matern52')
                            d_s2e{n_psi-n_time+p+1}=(-sqrt(5)*obj.DistMat_p/par.theta_z^2-10/3.*obj.DistMat_p.^2/par.theta_z^3).*exp(-sqrt(5)*obj.DistMat_p/par.theta_z)+sigma_eta.*sqrt(5).*obj.DistMat_z/(par.theta_z^2);
                        else
                            d_s2e{n_psi-n_time+p+1}=sigma_eta.*obj.DistMat_p{1}/(par.theta_z(1)^2);
                            d_s2e{n_psi-n_time+p+2}=sigma_eta.*obj.DistMat_p{2}/(par.theta_z(2)^2);
                        end
                    else
                        blocks_base=[0 cumsum(dim(1)*ones(1,p))];
                        for b=1:length(blocks_base)-1
                            idx_block=blocks_base(b)+1:blocks_base(b+1);
                            if strcmp(par.correlation_type,'exponential')
                                d_s2e{n_psi-n_time+p+b}=sparse(rr,rr);
                                d_s2e{n_psi-n_time+p+b}(idx_block,idx_block)=sigma_eta(idx_block,idx_block).*obj.DistMat_z(idx_block,idx_block)/(par.theta_z(b)^2);
                            elseif strcmp(par.correlation_type,'matern32')
                                d_s2e{n_psi-n_time+p+b}=sparse(rr,rr);
                                d_s2e{n_psi-n_time+p+b}(idx_block,idx_block)=-sqrt(3)*obj.DistMat_z(idx_block,idx_block)/(par.theta_z(b)^2).*exp(-sqrt(3)*obj.DistMat_z(idx_block,idx_block)/par.theta_z(b))+sigma_eta(idx_block,idx_block).*sqrt(3).*obj.DistMat_z(idx_block,idx_block)/(par.theta_z(b)^2);
                            elseif strcmp(par.correlation_type,'matern52')
                                d_s2e{n_psi-n_time+p+b}=sparse(rr,rr);
                                d_s2e{n_psi-n_time+p+b}(idx_block,idx_block)=(-sqrt(5)*obj.DistMat_z(idx_block,idx_block)/par.theta_z(b)^2-10/3.*obj.DistMat_z(idx_block,idx_block).^2/par.theta_z(b)^3).*exp(-sqrt(5)*obj.DistMat_z(idx_block,idx_block)/par.theta_z(b))+sigma_eta(idx_block,idx_block).*sqrt(5).*obj.DistMat_z(idx_block,idx_block)/(par.theta_z(b)^2);
                            else
                                d_s2e{n_psi-n_time+p+b*2-1}=sparse(rr,rr);
                                d_s2e{n_psi-n_time+p+b*2}=sparse(rr,rr);
                                d_s2e{n_psi-n_time+p+b*2-1}(idx_block,idx_block)=sigma_eta(idx_block,idx_block).*obj.DistMat_z{1}(idx_block,idx_block)/(par.theta_z(1,b)^2);
                                d_s2e{n_psi-n_time+p+b*2}(idx_block,idx_block)=sigma_eta(idx_block,idx_block).*obj.DistMat_z{2}(idx_block,idx_block)/(par.theta_z(2,b)^2);
                            end
                        end
                    end
                    
                    %v_z
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        z=1;
                        for j=1:p
                            for i=j:p
                                Id=[];
                                Jd=[];
                                elements=[];
                                
                                result1=sigma_eta(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                                result2 = result1/par.v_z(i,j);
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=cat(1,Id,idx_I+blocks(j));
                                Jd=cat(1,Jd,idx_J+blocks(i));
                                elements=cat(1,elements,result2(L));
                                if i~=j
                                    Id=cat(1,Id,idx_J+blocks(i));
                                    Jd=cat(1,Jd,idx_I+blocks(j));
                                    elements=cat(1,elements,result2(L));
                                end
                                
                                zero_density=(1-length(Id)/(N^2))*100;
                                if obj.tapering||zero_density>90
                                    d_s2e{n_psi-n_time+p+n_z_theta+z}=sparse(Id,Jd,elements,Np,Np);
                                else
                                    d_s2e{n_psi-n_time+p+n_z_theta+z}=zeros(Np);
                                    temp=sub2ind([Np,Np],Id,Jd);
                                    d_s2e{n_psi-n_time+p+n_z_theta+z}(temp)=elements;
                                end
                                z=z+1;
                            end
                        end
                    end
                    
                    if obj.stem_data.stem_modeltype.is('f-HDGM')
                        z=1;
                        for j=1:p
                            Id=[];
                            Jd=[];
                            elements=[];
                            
                            result1=sigma_eta(blocks(j)+1:blocks(j+1),blocks(j)+1:blocks(j+1));
                            result2 = result1/par.v_z(j,j);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=cat(1,Id,idx_I+blocks(j));
                            Jd=cat(1,Jd,idx_J+blocks(j));
                            elements=cat(1,elements,result2(L));
                           
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>90
                                d_s2e{n_psi-n_time+p+n_z_theta+z}=sparse(Id,Jd,elements,Np,Np);
                            else
                                d_s2e{n_psi-n_time+p+n_z_theta+z}=zeros(Np);
                                temp=sub2ind([Np,Np],Id,Jd);
                                d_s2e{n_psi-n_time+p+n_z_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end

                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  preliminary computations for time-invariant case  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(data.X_tv)&&size(par.sigma_eps,3)==1
                for i=1:n_psi 
                    if par.flag_sigma_eps_spline==1
                        d_Sgeo{i}=sparse(n_eps,1);
                    else
                        d_Sgeo{i}=sparse(N,N);
                    end
                end

                if not(isempty(data.X_bp))
                    a_bp_X=stem_misc.D_apply(data.X_bp{1},aj_bp,'l');
                    for i=1:n_psi
                        temp=zeros(N,1);
                        for j=1:n_bp_alpha
                            temp=temp+stem_misc.D_apply(data.X_bp{1},d_alpha_bp{i,j},'l');
                        end
                        d_a_bp_X=sparse(temp);
                        temp=stem_misc.M_apply(sigma_W_b,M,'b');
                        temp=stem_misc.D_apply(temp,d_a_bp_X,'l');
                        temp=stem_misc.D_apply(temp,a_bp_X,'r');
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                        
                        temp=stem_misc.D_apply(d_M_Sb{i},a_bp_X,'b');
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                        
                        temp=stem_misc.M_apply(sigma_W_b,M,'b');
                        temp=stem_misc.D_apply(temp,d_a_bp_X,'r');
                        temp=stem_misc.D_apply(temp,a_bp_X,'l');
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                    end
                end
                if not(isempty(data.X_p))
                    for z=1:k
                        aj_p=[ones(Np,1);zeros(N-Np,1)];
                        a_p_X=stem_misc.D_apply(data.X_p{1}(:,z),aj_p,'l');
                        for i=1:n_psi
                            temp=stem_misc.D_apply(d_Sp{i,z},a_p_X,'b');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,90);
                        end
                    end
                end
                   
                if sigma_eps_estimated
                    for i=1:n_psi
                        d_Sgeo{i}=d_Sgeo{i}+d_Seps{i};
                    end
                end
                
                if not(isempty(data.X_z))
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        temp=data.X_z{1};
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        if obj.stem_data.stem_modeltype.is('f-HDGM')
                            X_z_orlated=[data.X_z{1};zeros(N-size(data.X_z{1},1),size(data.X_z{1},2))];
                        else
                            X_z_orlated=[data.X_z{1};zeros(N-size(data.X_z{1},1),size(data.X_z{1},2))];
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Information matrix evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp('Hessian evaluation...');
            c0 = 0;
            IM=zeros(n_psi);
            if exit_toll>0
                last_varcov_t=zeros(n_psi);
            end
            tot=n_psi*(n_psi+1)/2*T;
            counter=1;
            
            for t=2:T+1
                if data.X_bp_tv
                    tBP=t-1;
                else
                    tBP=1;
                end
                if data.X_z_tv
                    tT=t-1;
                else
                    tT=1;
                end
                if data.X_beta_tv
                    tbeta=t-1;
                else
                    tbeta=1;
                end
                if data.X_p_tv
                    tP=t-1;
                else
                    tP=1;
                end
                if size(par.sigma_eps,3)>1
                    tS=t-1;
                else
                    tS=1;
                end
                
                if par.flag_sigma_eps_spline==1
                    tSeps=t-1;
                else
                    tSeps=1;
                end
                
                if data.X_tv||size(par.sigma_eps,3)>1
                    %compute sigma_geo in the time-variant case
                    if not(isempty(data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),data.X_bp{tBP},'b'),aj_bp,'b');
                    end
                    if not(isempty(data.X_p))
                        if isempty(data.X_bp)
                            if (obj.tapering)
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for z=1:k
                            aj_p=[ones(Np,1);zeros(N-Np,1)];
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{z},data.X_p{tP}(:,z),'b'),aj_p,'b');
                        end
                    end
                    
                    if isempty(data.X_p)&&isempty(data.X_bp)
                        sigma_geo=sigma_eps(:,:,tS);
                    else
                        sigma_geo=sigma_geo+sigma_eps(:,:,tS);
                    end
                    
                    if not(isempty(data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),data.X_bp{tBP},'b'),aj_bp,'b');
                    end
                    if not(isempty(data.X_p))
                        if not(isempty(sigma_geo))
                            if obj.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for z=1:k
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{z},data.X_p{tP}(:,z),'b'),aj_p,'b');
                        end
                    end
                   
                    if not(isempty(sigma_geo))
                        sigma_geo=sigma_eps(:,:,tS);
                    else
                        sigma_geo=sigma_geo+sigma_eps(:,:,tS);
                    end
                    
                    %compute d_Sgeo in the time-variant case
                    for i=1:n_psi 
                        if par.flag_sigma_eps_spline==1
                            d_Sgeo{i}=sparse(n_eps,1);
                        else
                            d_Sgeo{i}=sparse(N,N);
                        end 
                    end
                    if not(isempty(data.X_bp))
                        a_bp_X=data.X_bp{tBP};
                        for i=1:n_psi
                            temp=stem_misc.D_apply(d_M_Sb{i},a_bp_X,'b');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,90);
                        end
                    end
                    if not(isempty(data.X_p))
                        for z=1:k
                            a_p_X=data.X_p{tP}(:,z);
                            for i=1:n_psi
                                temp=stem_misc.D_apply(d_Sp{i,z},a_p_X,'b');
                                d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,90);
                            end
                        end
                    end
                    if sigma_eps_estimated
                        for i=1:n_psi
                            d_Sgeo{i}=d_Sgeo{i}+d_Seps{i,tS};
                        end
                    end
                    if not(isempty(data.X_z))
                        if obj.stem_data.stem_modeltype.is('HDGM')
                            temp=data.X_z{tT};
                            temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[data.X_z{tT};zeros(N-size(data.X_z{tT},1),size(data.X_z{tT},2))];
                        end
                    end
                end
                
                Lt=not(isnan(data.Y(:,t-1)));
                if t>2
                    Lt_lag=not(isnan(data.Y(:,t-2)));
                end
                if not(isempty(data.X_z))
                    if obj.stem_data.stem_modeltype.is('HDGM')
                        temp=data.X_z{tT};
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        if obj.stem_data.stem_modeltype.is('f-HDGM')
                            X_z_orlated=[data.X_z{tT};zeros(N-size(data.X_z{tT},1),size(data.X_z{tT},2))];
                        else
                            X_z_orlated=[data.X_z{tT};zeros(N-size(data.X_z{tT},1),size(data.X_z{tT},2))];
                        end
                    end
                    
                    if n_beta>0
                        X_beta_orlated=data.X_beta{tbeta};
                        X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                        e_t_Lt=data.Y(Lt,t-1)-X_beta_orlated(Lt,:)*par.beta-X_z_orlated(Lt,:)*zk_f(:,t);
                    else
                        e_t_Lt=data.Y(Lt,t-1)-X_z_orlated(Lt,:)*zk_f(:,t);
                    end
                    
                    if size(sigma_geo,3)==1
                        if iscell(sigma_geo) 
                            sigma_t_Lt=X_z_orlated(Lt,:)*Pk_f{t}*X_z_orlated(Lt,:)'+sigma_geo{t-1}(Lt,Lt);
                        else
                            sigma_t_Lt=X_z_orlated(Lt,:)*Pk_f{t}*X_z_orlated(Lt,:)'+sigma_geo(Lt,Lt);
                        end  
                    else
                        sigma_t_Lt=X_z_orlated(Lt,:)*Pk_f{t}*X_z_orlated(Lt,:)'+sigma_geo(Lt,Lt,t-1);
                    end
                    
                    d_Sz=cell(n_psi,1);
                    for i=1:n_psi
                        d_Sz{i}=sparse(N,N);
                    end
                else
                    if n_beta>0
                        X_beta_orlated=data.X_beta{tbeta};
                        X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                        e_t_Lt=data.Y(Lt,t-1)-X_beta_orlated(Lt,:)*par.beta;
                    else
                        e_t_Lt=data.Y(Lt,t-1);
                    end
                    if size(sigma_geo,3)==1
                        sigma_t_Lt=sigma_geo(Lt,Lt);
                    else
                        sigma_t_Lt=sigma_geo(Lt,Lt,t-1);
                    end
                end
                
                sigma_t_Lt=stem_misc.sparseif(sigma_t_Lt,90);

                for i=1:n_psi
                    if t==2
                        if not(isempty(data.X_z))
                            d_P{i}=d_s2e{i};

                            temp=X_z_orlated*d_P{i}*X_z_orlated';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,90);

                            
                            if par.flag_sigma_eps_spline==1
                                basis_Seps = full(getbasismatrix(data.X_h(:,tSeps),data.stem_fda.spline_basis_sigma));
                                if par.flag_logsigma==1
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+ sigma_geo{t-1}(Lt,Lt)*diag(basis_Seps(Lt,:)*d_Sgeo{i});
                                elseif par.flag_sqrsigma==1
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+ 2.*sqrt(sigma_geo{t-1}(Lt,Lt))*diag(basis_Seps(Lt,:)*d_Sgeo{i}); 
                                else
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+diag(basis_Seps(Lt,:)*d_Sgeo{i});
                                end
                            else 
                                d_St_Lt{i}=d_Sz{i}(Lt,Lt)+d_Sgeo{i}(Lt,Lt);
                            end
                            
                            d_J_Lt{i}=(d_G{i}*Pk_f{t}*X_z_orlated(Lt,:)'+G*d_P{i}*X_z_orlated(Lt,:)'-J{t}(:,Lt)*d_St_Lt{i})/sigma_t_Lt;
                            d_Z{i}=d_G{i}*zk_f(:,t-1);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        end
                    else
                        if not(isempty(data.X_z))
                            d_P{i}=d_G{i}*Pk_f{t-1}*G+G*d_P_lag{i}*G'+G*Pk_f{t-1}*d_G{i}'+d_s2e{i}-d_J_Lt_lag{i}*sigma_t_Lt_lag*J{t-1}(:,Lt_lag)'-J{t-1}(:,Lt_lag)*d_St_Lt_lag{i}*J{t-1}(:,Lt_lag)'-J{t-1}(:,Lt_lag)*sigma_t_Lt_lag*d_J_Lt_lag{i}';
                            
                            temp=X_z_orlated*d_P{i}*X_z_orlated';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,90);

                            if par.flag_sigma_eps_spline==1
                                basis_Seps = full(getbasismatrix(data.X_h(:,tSeps),data.stem_fda.spline_basis_sigma));
                                if par.flag_logsigma==1
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+ sigma_geo{t-1}(Lt,Lt)*diag(basis_Seps(Lt,:)*d_Sgeo{i});
                                elseif par.flag_sqrsigma==1
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+ 2.*sqrt(sigma_geo{t-1}(Lt,Lt))*diag(basis_Seps(Lt,:)*d_Sgeo{i}); 
                                else
                                    d_St_Lt{i}=d_Sz{i}(Lt,Lt)+diag(basis_Seps(Lt,:)*d_Sgeo{i});
                                end
                            else
                                d_St_Lt{i}=d_Sz{i}(Lt,Lt)+d_Sgeo{i}(Lt,Lt);
                            end
                            
                            d_J_Lt{i}=(d_G{i}*Pk_f{t}*X_z_orlated(Lt,:)'+G*d_P{i}*X_z_orlated(Lt,:)'-J{t}(:,Lt)*d_St_Lt{i})/sigma_t_Lt;
                            d_Z{i}=d_G{i}*zk_f(:,t-1)+G*d_Z_lag{i}+d_J_Lt_lag{i}*e_t_Lt_lag+J{t-1}(:,Lt_lag)*d_e_Lt_lag{i};
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i}-X_z_orlated(Lt,:)*d_Z{i};
                            else
                                d_e_Lt{i}=-X_z_orlated(Lt,:)*d_Z{i};
                            end
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        end
                    end
                end
                
                temp0=cell(n_psi,1);
                temp1=cell(n_psi,1);
                
                for i=1:n_psi
                    if issparse(sigma_t_Lt)
                        r = symamd(sigma_t_Lt);
                        c=chol(sigma_t_Lt(r,r));
                        temp0{i}(r,:)=stem_misc.chol_solve(c,d_e_Lt{i}(r,:));
                    else
                        c=chol(sigma_t_Lt);
                        temp0{i}=stem_misc.chol_solve(c,d_e_Lt{i});
                    end
                    
                    d_St_i_Lt=d_St_Lt{i};
                    if nnz(d_St_i_Lt)>0
                        if issparse(sigma_t_Lt)
                            temp1{i}(r,:)=full(stem_misc.chol_solve(c,d_St_i_Lt(r,:)));
                        else
                            temp1{i}=full(stem_misc.chol_solve(c,d_St_i_Lt));
                        end
                    else
                        temp1{i}=spalloc(size(sigma_t_Lt,1),size(sigma_t_Lt,2),0);
                    end
                end
                
                for i=1:n_psi
                    blocks=0:50:size(temp1{i},1);
                    if not(blocks(end)==size(temp1{i},1))
                        blocks=cat(2,blocks,size(temp1{i},1));
                    end
                    nnz_temp1_i=nnz(temp1{i});
                    
                    for j=i:n_psi
                        IM(i,j)=IM(i,j)+d_e_Lt{i}'*temp0{j};
                        if (nnz_temp1_i>0)&&(nnz(temp1{j})>0)
                            sumtrace=0;
                            for z=1:length(blocks)-1
                                idx=blocks(z)+1:blocks(z+1);
                                sumtrace=sumtrace+trace(temp1{i}(idx,:)*temp1{j}(:,idx));
                            end
                            IM(i,j)=IM(i,j)+0.5*sumtrace+0.25*trace(temp1{i})*trace(temp1{j});
                        end
                        counter=counter+1;
                        if (mod(counter,100)==0)||(counter==tot)
                            if mod(round(counter/tot*100),1)==0&&c0~=round(counter/tot*100)
                                c0 = round(counter/tot*100);
                                disp(['Hessian evaluation: ',num2str(round(counter/tot*100)),'% completed']);
                            end
                        end
                    end
                end
                
                if not(isempty(data.X_z))
                    d_P_lag=d_P;
                    d_J_Lt_lag=d_J_Lt;
                    d_Z_lag=d_Z;
                end
                
                d_e_Lt_lag=d_e_Lt;
                e_t_Lt_lag=e_t_Lt;
                d_St_Lt_lag=d_St_Lt;
                sigma_t_Lt_lag=sigma_t_Lt;
                if data.X_tv
                    sigma_geo=[];
                end
                
                if (exit_toll>0)&&(t>2)
                    temp_IM=IM+triu(IM,1)';
                    current_varcov_t=inv(temp_IM/t*T);
                    delta=norm(current_varcov_t-last_varcov_t,'fro')/norm(current_varcov_t,'fro');
                    if delta<exit_toll
                        IM=IM+triu(IM,1)';
                        obj.stem_par.varcov=inv(IM/t*T);
                        obj.stem_EM_result.stem_par.varcov=obj.stem_par.varcov;
                        approximated=1;
                        disp('Exit tolerance reached. Aproximated VarCov computed.');
                        break
                    else
                        approximated=0;
                    end
                    last_varcov_t=current_varcov_t;
                else
                    approximated=0;
                end
            end
            if approximated==0
                IM=IM+triu(IM,1)';
                obj.stem_par.varcov=inv(IM);
                obj.stem_EM_result.stem_par.varcov=obj.stem_par.varcov;
            end
        end   
        
        function set_initial_values(obj,stem_par)
            %DESCRIPTION: set the initial values of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %none: the stem_par_initial property is updated                
            if not(isa(stem_par,'stem_par'))
                error('The input argument must be of class stem_par');
            end
            obj.stem_par_initial=stem_par;
            if obj.stem_par_initial.p>0
                if obj.stem_par.stem_par_constraints.time_diagonal&&not(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'}))
                    obj.stem_par_initial.G=diag(diag(obj.stem_par_initial.G));
                    obj.stem_par_initial.sigma_eta=diag(diag(obj.stem_par_initial.sigma_eta));
                end
            end
            if not(isempty(stem_par.G))
                G=stem_par.G;
                temp=abs(G-eye(size(G)));
                if sum(temp(:))==0
                    warning('Matrix G is the identity matrix. The temporal component z is a random walk and G will not be updated by the EM algorhithm');
                end
            end
            
            if not(isempty(stem_par.v_p))
                if size(stem_par.v_p,3)==1&&stem_par.k>1
                    warning(['The same cross-correlation matrix v_p is used to modell all the ',num2str(stem_par.k),' latent spatial variables w_p']);
                end
            end
            obj.stem_par=obj.stem_par_initial;
        end
         
        function set_distance_matrix(obj,type,force)
            %DESCRIPTION: generates the distance matrices
            %
            %INPUT
            %obj          - [stem_data object] (1x1) the stem_data object
            %<force>      - [integer]          (1x1) (Default: 0) 1: force to compute the distance matrix even if the model does not include spatial latent variables 
            %<type>       - [string]           (1x1) (Default: 'both') 'point': only the distance matrix for the point data is evaluated. 
            %                                                          'pixel': only the distance matrix for the pixel data is evaluated.
            %                                                          'both':  both the matrices are evaluated.
            %OUTPUT         
            %
            %none: the DistMat_p and DistMat_b property are generated and updated

            if nargin<2
                type='both';
            end
            if nargin<3
                force=0;
            end

            cmp=strcmp(type,{'both','point','pixel'});
            if sum(cmp)==0
                error('type must be point, pixel or both');
            end
            
            if strcmp(type,'point')||strcmp(type,'both')
                if not(isempty(obj.stem_data.stem_varset_p.X_p))||force||(obj.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})&&not(isempty(obj.stem_data.stem_varset_p.X_z)))
                    disp('Generating point distance matrices...');
                    if not(obj.stem_data.stem_modeltype.is('f-HDGM'))
                        obj.DistMat_p=obj.stem_data.stem_gridlist_p.get_distance_matrix(obj.stem_data.stem_modeltype.model_name,obj.stem_par.correlation_type,1,[]);
                    else
                        if obj.stem_data.stem_modeltype.is('f-HDGM')
                            %y_t and the latent variable z_t have different dimensions
                            type=1;
                            idx_var=1;
                            DistMat=obj.stem_data.stem_gridlist_p.get_distance_matrix(obj.stem_data.stem_modeltype.model_name,obj.stem_par.correlation_type,type,idx_var);
                            if iscell(DistMat)
                                for i=1:length(DistMat)
                                    obj.DistMat_z{i}=kron(ones(size(obj.stem_data.stem_varset_p.X_z{1},2)),DistMat{i});
                                end
                            else
                                obj.DistMat_z=kron(ones(size(obj.stem_data.stem_varset_p.X_z{1},2)),DistMat);
                            end
                        end
                    end
                    disp('Generation ended.');
                end
            end
            if strcmp(type,'pixel')||strcmp(type,'both')
                if not(isempty(obj.stem_data.stem_gridlist_b))&&not(isempty(obj.stem_data.stem_varset_b.X_bp))
                    disp('Generating pixel data distance matrices...');
                    obj.DistMat_b=obj.stem_data.stem_gridlist_b.get_distance_matrix(obj.stem_data.stem_modeltype.model_name,obj.stem_par.correlation_type,obj.stem_par.stem_par_constraints.pixel_correlated,[]);
                    disp('Generation ended.');
                end
            end
        end     

        %Initial values estimation functions 
        function [beta0] = get_beta0(obj)
            %DESCRIPTION: provides the initial values of the beta parameter
            %
            %INPUT
            %
            %obj               - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %
            %beta0             - [double] (kx1) the initial beta parameter
            if obj.stem_par.n_beta>0
                N = size(obj.stem_data.X_beta{1},1);
                y = obj.stem_data.Y(1:N,:);
                y=y(:);
                T = obj.T;
                x = zeros(N*T,size(obj.stem_data.X_beta{1},2));
                
                for t=1:T
                    if length(obj.stem_data.X_beta)==T
                        tT=t;
                    else
                        tT=1;
                    end
                    x((t-1)*N+1:t*N,:)=obj.stem_data.X_beta{tT};
                end
                
                L=not(isnan(y));
                beta0 = (x(L,:)'*x(L,:))\x(L,:)'*y(L);
            else
                disp('WARNING: the model does not include data to estimate beta');
                beta0=[];
            end
        end
        
        function coe_sigma_eps0 = get_coe_log_sigma_eps0(obj)
            %DESCRIPTION: provides the initial values of the coefficients of basis to estimate sigma_eps2(h)
            %
            %INPUT
            %
            %obj             - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT       
            %
            %coe_sigma_eps0  - [double] (nbasis_sigmax1) the initial coefficients of basis when estimating sigma_eps2(h)
            if obj.stem_par.flag_sigma_eps_spline~=1
                error('The function is only available in the case that sigma_eps2(h) is estimated for the f-HDGM ');
            end
            T = obj.T;
            N = size(obj.stem_data.X_h,1);
            res = zeros(N,T);
            beta0 = obj.get_beta0();
            for t=1:T 
                if length(obj.stem_data.X_beta)==T
                    tT=t;
                else
                    tT=1;
                end
                res(1:N,t) = obj.stem_data.Y(:,t)-obj.stem_data.X_beta{tT}*beta0;
            end

            if obj.stem_par.flag_logsigma==1
                y = res(1:N,:).^2;
                X_h=obj.stem_data.X_h(:);
                L1=not(isnan(X_h));
                bmat = full(getbasismatrix(X_h,obj.stem_data.stem_fda.spline_basis_sigma));
                y=y(:);
                L2=not(isnan(y));
                L=L1&L2;
                initialvalue = (bmat(L,:)'*bmat(L,:))\bmat(L,:)'*log(y(L));
                coe_sigma_eps0 = initialvalue;

            else
                y=res(1:N,:).^2;
                y = y(:);
                T = obj.T;
                x = zeros(N*T,obj.stem_par.k_sigma);
                for t=1:T
                    if length(obj.stem_data.X_h)==T
                        tT=t;
                    else
                        tT=1;
                    end
                    x((t-1)*N+1:t*N,:)=full(getbasismatrix(obj.stem_data.X_h(:,tT),obj.stem_data.stem_fda.spline_basis_sigma));
                end
                L=not(isnan(y));
                coe_sigma_eps0 = (x(L,:)'*x(L,:))\x(L,:)'*y(L);
            end
        end
        
        %Export functions. Useful to avoid access to the properties of nested objects
        function N = N(obj)
            N=obj.stem_data.N();
        end
        
        function Nb = Nb(obj)
            Nb=obj.stem_data.Nb();
        end
        
        function Np = Np(obj)
            Np=obj.stem_data.Np();
        end
        
        function T = T(obj)
            T=obj.stem_data.T();
        end
        
        function nvar=nvar(obj)
            nvar=obj.stem_data.nvar();
        end
        
        function dim=dim(obj)
            dim=obj.stem_data.dim();
        end
        
        %Class set functions 
        function set.stem_data(obj,stem_data)
            if isa(stem_data,'stem_data')
                obj.stem_data=stem_data;
            else
                error('The argument must be of class stem_data');
            end
        end
        
        function set.stem_par(obj,stem_par)
            if not(isa(stem_par,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            end
            obj.stem_par=stem_par;
        end        
       
        function set.stem_par_initial(obj,stem_par_initial)
            if not(isa(stem_par_initial,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            end
            obj.stem_par_initial=stem_par_initial;
        end
    end
end
