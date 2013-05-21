%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-STEM - Distributed Space Time Expecation Maximization      %
%                                                              %
% Author: Francesco Finazzi                                    %
% E-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo - Dept. of Engineering    %
% Author website: http://www.unibg.it/pers/?francesco.finazzi  %
% Code website: https://code.google.com/p/d-stem/              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_model < handle
    
    %CONSTANTS
    %N  = n1_g+...+nq_g+n1_r+...+nq_r - total number of observation sites
    %N_g = n1_g+...+nq_g - total number of point sites
    %N_r = n1_r+...+nq_r - total number of pixel sites
    %S = 1 if only point data are considered and S=2 if both point and pixel data are considered
    
    properties
        stem_data=[];           %[stem_data object] (1x1) object containing all the data used to estimated the model
        stem_par=[];            %[stem_par object]  (1x1) parameter set updated at each iteration of the EM algorithm
        stem_par_initial=[];    %[stem_par object]  (1x1) starting parameter set
        stem_par_sim=[];        %[stem_par object]  (1x1) parameter set used to simulate data (if data are simulated)
        estimated=0;            %[boolean] (1x1) 0: the model has not been estimated; 1: the model has been estimated
    end
    
    properties (SetAccess = private)
        stem_EM_result=[];      %[stem_EM_result object] (1x1) object containing all the results of the EM estimation
        cross_validation=0;     %[boolean] (1x1) 0: the model has been estimated considering all the data; 1: the model has bee estimated excluding the cross-validation data.
        system_size=100;        %[integer] (1x1) if N > system_size than only the diagonal is computed in matrix multiply operations
        tapering=[];            %[boolean] (1x1) 0:tapering is not enabled; 1:tapering is enabled on point sites or pixel sites
        tapering_r=[];          %[boolean] (1x1) 0:tapering is not enabled on pixel sites; 1:tapering is enabled on pixel sites
        tapering_g=[];          %[boolean] (1x1) 0:tapering is not enabled on point sites; 1:tapering is enabled on point sites
    end
    
    methods
        
        function obj = stem_model(stem_data,stem_par)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_data - [stem_data object]  (1x1)
            %stem_par  - [stem_par object]   (1x1)
            %
            %OUTPUT
            %obj       - [stem_model object] (1x1)
            if nargin<2
                error('Not enough input arguments');
            end
            obj.stem_data=stem_data;
            obj.stem_par=stem_par;
            
            if not(isempty(stem_data.stem_gridlist_r))
                if not(isempty(stem_data.stem_gridlist_r.tap))
                    obj.tapering_r=1;
                else
                    obj.tapering_r=0;
                end
            end
            if not(isempty(stem_data.stem_gridlist_g.tap))
                obj.tapering_g=1;
            else
                obj.tapering_g=0;
            end
            if not(isempty(stem_data.stem_gridlist_r))
                obj.tapering=obj.tapering_g|obj.tapering_r;
            else
                obj.tapering=obj.tapering_g;
            end
            
            if stem_par.theta_clustering>0
                stem_data.update_distance('both',1);
            end
        end
        
        function [aj_rg,aj_g] = get_aj(obj)
            %DESCRIPTION: provides the vector aj_rg and aj_g used in the EM estimation
            %
            %INPUT
            %obj   - [stem_model object] (1x1)
            %
            %OUTPUT
            %aj_rg - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_rg
            %aj_g  - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_g. 
            
            %NOTE
            %The elements of aj_g from Ng+1 to N are all zeros. This allows the
            %use of stem_misc.D_apply both for the pixel data and the point
            %level data avoiding the use of J_rg and J_g
            if not(isempty(obj.stem_data.stem_varset_r))
                aj_rg=zeros(obj.stem_data.N,1);
                blocks=[0 cumsum(obj.stem_data.dim)];
                for i=1:obj.stem_data.nvar
                    aj_rg(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_rg(i);
                end
            else
                aj_rg=[];
            end
            if obj.stem_par.k>0
                aj_g=zeros(obj.stem_data.N,obj.stem_par.k);
                blocks=[0 cumsum(obj.stem_data.stem_varset_g.dim)];
                for k=1:obj.stem_par.k
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        aj_g(blocks(i)+1:blocks(i+1),k)=obj.stem_par.alpha_g(i,k);    
                    end
                end
            else
                aj_g=[];
            end
        end
        
        function [aj_rg_r,j_r] = get_jrg(obj,r)
            %DESCRIPTION: provides the vectors aj_rg_r and j_r used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_rg_r - [double]            (Nx1) is the vector with elements equal to alpha_rg(r) only for the sites of the r-th variable
            %j_r     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            aj_rg_r=zeros(obj.stem_data.N,1);
            j_r=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.nvar
                if i==r
                    j_r(blocks(i)+1:blocks(i+1))=1;
                    aj_rg_r(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_rg(i);
                end
            end
        end
        
        function [aj_g_rs,j_r] = get_jg(obj,r,s)
            %DESCRIPTION: provides the vectors aj_g_rs and j_r used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %s       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_g_rs - [double]            (Nx1) is the vector with elements equal to alpha_g(i,s) only for the sites of the r-th variable
            %j_r     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            aj_g_rs=zeros(obj.stem_data.N,1);
            j_r=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.stem_varset_g.nvar
                if i==r
                    j_r(blocks(i)+1:blocks(i+1))=1;
                    aj_g_rs(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_g(i,s); 
                end
            end
        end        
        
        function [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = get_sigma(obj,sigma_W_r)
            %DESCRIPTION: provides the variance-covariance matrices and some vectors that are used in the EM algorithm
            %
            %INPUT
            %
            %obj         - [stem_model object] (1x1)
            %<sigma_W_r> - [double]            (NrxNr) (default: []) variance-covariance matrix of W_r. It is provided as input argument during kriging when sigma_W_r does not change across blocks  
            %
            %OUTPUT
            %
            %sigma_eps   - [double]            (NxN) variance-covariance matrix of epsilon
            %sigma_W_r   - [double]            (N_rxN_r) variance-covariance matrix of W_r
            %sigma_W_g   - [double]            {K}(N_gxN_g) variance-covariance matrices of the K W_g_i
            %sigma_geo   - [double]            (NxN) variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %sigma_Z     - [double]            (pxp) variance-covariance of Z
            %aj_rg       - [double]            (Nx1) see the details of the method get_aj;
            %aj_g        - [double]            (Nx1) see the details of the method get_aj;
            %M           - [integer >0]        (N_gx1) see the details of the method update_M of the class stem_data
            %
            %NOTE
            %sigma_geo is provided only if it is time-invariant otherwise it is evaluated at each step of the EM algorithm
            disp('    Marginal variance-covariance matrices evaluation started...');
            ct1=clock;
            if nargin==1
                sigma_W_r=[];
            end
            
            nvar=obj.stem_data.nvar;
            N=obj.stem_data.N;
            M=obj.stem_data.M;
            dim=obj.stem_data.dim;
            
            %sigma_eps
            d=[];
            for i=1:nvar
                d=[d;repmat(obj.stem_par.sigma_eps(i,i),dim(i),1)];
            end
            
            I=1:length(d);
            sigma_eps=sparse(I,I,d);
            
            %sigma_W_r
            if not(isempty(obj.stem_data.stem_varset_r))
                if not(isempty(obj.stem_data.stem_varset_r.X_rg))
                    if (nargin==1)||isempty(sigma_W_r)
                        if not(obj.tapering_r)
                            sigma_W_r=zeros(obj.stem_data.stem_varset_r.N);
                        end
                        blocks=[0 cumsum(obj.stem_data.stem_varset_r.dim)];
                        if obj.stem_par.pixel_correlated
                            if obj.tapering_r
                                I=zeros(nnz(obj.stem_data.DistMat_r),1);
                                J=zeros(nnz(obj.stem_data.DistMat_r),1);
                                elements=zeros(nnz(obj.stem_data.DistMat_r),1);
                            end
                            idx=0;
                            for j=1:obj.stem_data.stem_varset_r.nvar
                                for i=j:obj.stem_data.stem_varset_r.nvar
                                    idx_r=blocks(i)+1:blocks(i+1);
                                    idx_c=blocks(j)+1:blocks(j+1);
                                    if not(obj.tapering_r)
                                        sigma_W_r(idx_r,idx_c)=obj.stem_par.v_r(i,j)*stem_misc.correlation_function(...
                                            obj.stem_par.theta_r,obj.stem_data.DistMat_r(idx_r,idx_c),obj.stem_par.correlation_type);
                                        if not(i==j)
                                            sigma_W_r(idx_c,idx_r)=sigma_W_r(idx_r,idx_c)';
                                        end
                                    else
                                        corr_result=stem_misc.correlation_function(obj.stem_par.theta_r,obj.stem_data.DistMat_r(idx_r,idx_c),obj.stem_par.correlation_type);
                                        weights=stem_misc.wendland(obj.stem_data.DistMat_r(idx_r,idx_c),obj.stem_data.stem_gridlist_r.tap);
                                        corr_result.correlation=obj.stem_par.v_r(i,j)*corr_result.correlation.*weights;
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
                            if obj.tapering_r
                                sigma_W_r=sparse(I,J,elements);
                            end
                        else
                            if obj.tapering_r
                                idx=0;
                                nonzeros=0;
                                for i=1:obj.stem_data.stem_varset_r.nvar
                                    nonzeros=nonzeros+nnz(obj.stem_data.DistMat_r(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)));
                                end
                                I=zeros(nonzeros,1);
                                elements=zeros(nonzeros,1);
                            end
                            for i=1:obj.stem_data.stem_varset_r.nvar
                                idx_rc=blocks(i)+1:blocks(i+1);
                                if not(obj.tapering_r)
                                    sigma_W_r(idx_rc,idx_rc)=stem_misc.correlation_function(obj.stem_par.theta_r(i,:),obj.stem_data.DistMat_r(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                else
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_r(i,:),obj.stem_data.DistMat_r(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                    weights=stem_misc.wendland(obj.stem_data.DistMat_r(idx_rc,idx_rc),obj.stem_data.stem_gridlist_r.tap);
                                    corr_result.correlation=obj.stem_par.v_r(i,i)*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                            if obj.tapering_r
                                sigma_W_r=sparse(I,J,elements);
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
            
           %sigma_W_g
           if obj.stem_par.k>0
               blocks=[0 cumsum(obj.stem_data.stem_varset_g.dim)];
               for k=1:obj.stem_par.k
                   if obj.tapering_g
                       I=zeros(nnz(obj.stem_data.DistMat_g),1);
                       J=zeros(nnz(obj.stem_data.DistMat_g),1);
                       elements=zeros(nnz(obj.stem_data.DistMat_g),1);
                   else
                       sigma_W_g{k}=obj.stem_data.DistMat_g;   
                   end
                   idx=0;
                   for j=1:obj.stem_data.stem_varset_g.nvar
                       for i=j:obj.stem_data.stem_varset_g.nvar
                           idx_r=blocks(i)+1:blocks(i+1);
                           idx_c=blocks(j)+1:blocks(j+1);
                           if not(obj.tapering_g)
                               sigma_W_g{k}(idx_r,idx_c)=obj.stem_par.v_g(i,j,k)*stem_misc.correlation_function(...
                                   obj.stem_par.theta_g(k),obj.stem_data.DistMat_g(idx_r,idx_c),obj.stem_par.correlation_type);
                               if not(i==j)
                                   sigma_W_g{k}(idx_c,idx_r)=sigma_W_g{k}(idx_r,idx_c)';
                               end
                           else
                               corr_result=stem_misc.correlation_function(obj.stem_par.theta_g(k),obj.stem_data.DistMat_g(idx_r,idx_c),obj.stem_par.correlation_type);
                               weights=stem_misc.wendland(obj.stem_data.DistMat_g(idx_r,idx_c),obj.stem_data.stem_gridlist_g.tap);
                               corr_result.correlation=obj.stem_par.v_g(i,j,k)*corr_result.correlation.*weights;
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
                   if obj.tapering_g
                       sigma_W_g{k}=sparse(I,J,elements);
                   end
               end
           else
               sigma_W_g=[];
           end
           clear I
           clear J
           clear elements
           clear weights
           clear corr_result
           
           %sigma_geo
           [aj_rg,aj_g]=obj.get_aj;
           if not(obj.stem_data.X_tv)
               %time invariant case
               if not(isempty(obj.stem_data.stem_varset_r))
                   sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),obj.stem_data.X_rg(:,1,1),'b'),aj_rg,'b');
               end
               if obj.stem_par.k>0
                   if isempty(obj.stem_data.stem_varset_r)
                       %se manca il remoto allora sigma_geo non � stata
                       %ancora allocata
                       if obj.tapering
                           sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                       else
                           sigma_geo=zeros(N);
                       end
                   end
                   for k=1:obj.stem_par.k
                       sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},obj.stem_data.X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                   end
               end
               if (obj.stem_par.k==0)&&(isempty(obj.stem_data.stem_varset_r))
                   sigma_geo=sigma_eps;
               else
                   sigma_geo=sigma_geo+sigma_eps;
               end
           else
               sigma_geo=[]; %it means sigma_geo is time-variant and has to be evaluated at each time step
           end
           
           if obj.stem_par.p>0
               sigma_Z=(eye(obj.stem_par.p^2)-kron(obj.stem_par.G,obj.stem_par.G))\obj.stem_par.sigma_eta(:);%variance of a VAR (Lutkepohl pag.27)
               sigma_Z=reshape(sigma_Z,obj.stem_par.p,obj.stem_par.p);
           else
               sigma_Z=[];
           end   
           ct2=clock;
           disp(['    Marginal variance-covariance matrices evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            %DESCRIPTION: is the front-end method to simulate data from this stem_model object
            %
            %INPUT
            %
            %obj               - [stem_model object] (1x1)
            %<nan_rate>        - [double [0,1]]      (Sx1) (default: []) missing rate of the simulated data
            %<nan_pattern_par> - [double >0]         (Sx1) (default: []) (UoM: km) parameter of the exponential spatial correlation function used to define the spatial pattern of the missing data
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
            %obj                - [stem_model object]      (1x1)
            %stem_EM_options    - [stem_EM_options object] (1x)
            %
            %OUTPUT
            %
            %none: the stem_EM object is created. The method stem_EM.estimate is used to estimate the model and it updates the stem_par object

            if not(isempty(obj.stem_data.stem_crossval))
                disp('Data modification for cross-validation started...');
                idx_var=obj.stem_data.stem_varset_g.get_Y_index(obj.stem_data.stem_crossval.variable_name);
                if isempty(idx_var)
                    error('Variable not found. Has it been deleted?');
                end
                
                %recover the indices of the cross-validation sites 
                indices=obj.stem_data.stem_crossval.indices;
                Y={obj.stem_data.stem_varset_g.Y{idx_var}(indices,:)};
                Y_name={obj.stem_data.stem_varset_g.Y_name{idx_var}};
                if not(isempty(obj.stem_data.stem_varset_g.X_rg))
                    X_rg={obj.stem_data.stem_varset_g.X_rg{idx_var}(indices,:,:)};
                    X_rg_name={obj.stem_data.stem_varset_g.X_rg{idx_var}};
                else
                    X_rg={};
                    X_rg_name={};
                end
                if not(isempty(obj.stem_data.stem_varset_g.X_beta))
                    X_beta={obj.stem_data.stem_varset_g.X_beta{idx_var}(indices,:,:)};
                    X_beta_name={obj.stem_data.stem_varset_g.X_beta_name{idx_var}};
                else
                    X_beta={};
                    X_beta_name={};
                end 
                if not(isempty(obj.stem_data.stem_varset_g.X_z))
                    X_z={obj.stem_data.stem_varset_g.X_z{idx_var}(indices,:,:)};
                    X_z_name={obj.stem_data.stem_varset_g.X_z_name{idx_var}};
                else
                    X_z={};
                    X_z_name={};
                end          
                if not(isempty(obj.stem_data.stem_varset_g.X_g))
                    X_g={obj.stem_data.stem_varset_g.X_g{idx_var}(indices,:,:,:)};
                    X_g_name={obj.stem_data.stem_varset_g.X_g_name{idx_var}};
                else
                    X_g={};
                    X_g_name={};
                end      
                
                %set the cross_mindistance vector
                if not(isempty(obj.stem_data.DistMat_g))
                    dim=obj.stem_data.dim;
                    blocks=[0 cumsum(dim)];
                    temp_dist=obj.stem_data.DistMat_g(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                    temp_dist=temp_dist(indices,:);
                    temp_dist(:,indices)=[];
                    obj.stem_data.stem_crossval.min_distance=min(temp_dist');
                    clear temp_dist
                end
                
                obj.stem_data.stem_crossval.stem_varset=stem_varset(Y,Y_name,X_rg,X_rg_name,X_beta,X_beta_name,X_z,X_z_name,X_g,X_g_name);
                obj.stem_data.stem_crossval.stem_gridlist=stem_gridlist();
                coordinate=obj.stem_data.stem_gridlist_g.grid{idx_var}.coordinate;
                st_grid=stem_grid(coordinate(indices,:),'deg','sparse','point');
                obj.stem_data.stem_crossval.stem_gridlist.add(st_grid);
                %remove the cross-validation data from the estimation dataset
                obj.stem_data.site_crop(obj.stem_data.stem_crossval.type,obj.stem_data.stem_crossval.variable_name,indices);
                obj.cross_validation=1;
                disp('Data modification ended.');
            else
                obj.cross_validation=0;
            end
            
            st_EM=stem_EM(obj,stem_EM_options);
            %set the current parameter value with the estimated initial value
            obj.stem_par=obj.stem_par_initial;
            if isempty(stem_EM_options.pathparallel)
                obj.stem_EM_result=st_EM.estimate();
            else
                obj.stem_EM_result=st_EM.estimate_parallel(stem_EM_options.pathparallel);
            end
            obj.estimated=1;
            if obj.cross_validation
                st_krig=stem_krig(obj);
                block_size=1000;
                back_transform=0;
                no_varcov=1;
                crossval=1;
                obj.stem_data.stem_crossval.stem_krig_result=st_krig.kriging(obj.stem_data.stem_crossval.variable_name,obj.stem_data.stem_crossval.stem_gridlist.grid{1},block_size,[],[],back_transform,no_varcov,crossval);
                obj.stem_data.stem_crossval.res=obj.stem_data.stem_crossval.stem_varset.Y{idx_var}-obj.stem_data.stem_crossval.stem_krig_result.y_hat;
                obj.stem_data.stem_crossval.mse=nanvar(obj.stem_data.stem_crossval.res');
                obj.stem_data.stem_crossval.relative_mse=obj.stem_data.stem_crossval.mse./nanvar(obj.stem_data.stem_crossval.stem_varset.Y{idx_var}');

                s=obj.stem_data.stem_varset_g.Y_stds{idx_var};
                m=obj.stem_data.stem_varset_g.Y_means{idx_var};
                if (obj.stem_data.stem_varset_g.standardized)&&not(obj.stem_data.stem_varset_g.log_transformed)
                    y_hat_back=obj.stem_data.stem_crossval.stem_krig_result.y_hat*s+m;
                    y=obj.stem_data.stem_crossval.stem_varset.Y{idx_var}*s+m;
                    obj.stem_data.stem_crossval.res_backtransformed=y-y_hat_back;
                    obj.stem_data.stem_crossval.y_back=y;
                    obj.stem_data.stem_crossval.y_hat_back=y_hat_back;
                end
                if (obj.stem_data.stem_varset_g.standardized)&&(obj.stem_data.stem_varset_g.log_transformed)
                    y_hat_back=obj.stem_data.stem_crossval.stem_krig_result.y_hat;
                    y_hat_back=exp(y_hat_back*s+m+(s^2)/2);
                    y=exp(obj.stem_data.stem_crossval.stem_varset.Y{idx_var}*s+m);
                    obj.stem_data.stem_crossval.res_backtransformed=y-y_hat_back;
                    obj.stem_data.stem_crossval.y_back=y;
                    obj.stem_data.stem_crossval.y_hat_back=y_hat_back;
                end
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
                    if not(isempty(obj.tapering_g))
                        disp(['  Point data tapering: ',num2str(obj.stem_data.stem_gridlist_g.tap),' km']);
                    else
                        disp('   Tapering is NOT enabled on point data');
                    end
                    if not(isempty(obj.tapering_r))
                        disp(['  Pixel data tapering: ',num2str(obj.stem_data.stem_gridlist_r.tap),' km']);
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
                if not(isempty(obj.stem_par.beta))
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        disp(['* Beta coefficients related to the point variable ',obj.stem_data.stem_varset_g.Y_name{i}]);
                        output=[];
                        output{1,1}='Loading coefficient';
                        output{1,2}='Value';
                        output{1,3}='Std';
                        for j=1:length(obj.stem_data.stem_varset_g.X_beta_name{i})
                            output{j+1,1}=obj.stem_data.stem_varset_g.X_beta_name{i}{j};
                            output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{j+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                    end
                    if not(isempty(obj.stem_data.stem_varset_r))
                        if not(isempty(obj.stem_data.stem_varset_r.X_beta))
                            for i=1:obj.stem_data.stem_varset_r.nvar
                                disp(['* Beta coefficients related to the pixel variable ',obj.stem_data.stem_varset_r.Y_name{i}]);
                                output=[];
                                output{1,1}='Loading coefficient';
                                output{1,2}='Value';
                                output{1,3}='Std';
                                for j=1:length(obj.stem_data.stem_varset_r.X_beta_name{i})
                                    output{j+1,1}=obj.stem_data.stem_varset_r.X_beta_name{i}{j};
                                    output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                                    else
                                        output{j+1,3}='Not computed';
                                    end
                                    counter=counter+1;
                                end
                                disp(output);
                            end
                        end
                    end
                end
                output=[];
                disp('* Sigma_eps diagonal elements (Variance)')
                output{1,1}='Variable';
                output{1,2}='Value';
                output{1,3}='Std';
                for i=1:obj.stem_data.stem_varset_g.nvar
                    output{i+1,1}=obj.stem_data.stem_varset_g.Y_name{i};
                    output{i+1,2}=num2str(obj.stem_par.sigma_eps(i,i),'%05.3f');
                    if not(isempty(obj.stem_EM_result.varcov))
                        output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                    else
                        output{i+1,3}='Not computed';
                    end
                    counter=counter+1;
                end
                if not(isempty(obj.stem_data.stem_varset_r))
                    delta=obj.stem_data.stem_varset_g.nvar;
                    for i=1:obj.stem_data.stem_varset_r.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_r.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.sigma_eps(i+delta,i+delta),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                end
                disp(output);
                output=[];
                if not(isempty(obj.stem_data.stem_varset_r))
                    disp('* alpha_rg elements')
                    output{1,1}='Variable';
                    output{1,2}='Value';
                    output{1,3}='Std';
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        output{i+1,1}=obj.stem_data.stem_varset_g.Y_name{i};
                        output{i+1,2}=num2str(obj.stem_par.alpha_rg(i),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    delta=obj.stem_data.stem_varset_g.nvar;
                    for i=1:obj.stem_data.stem_varset_r.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_r.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.alpha_rg(i+delta),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    output=[];
                    if obj.stem_par.pixel_correlated
                        disp('* Pixel data are cross-correlated.');
                        disp(' ');
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        output{2,1}='Theta_r';
                        output{2,2}=num2str(obj.stem_par.theta_r(1),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{2,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{2,3}='Not computed';
                        end
                        counter=counter+1;
                        disp(output);
                        output=[];
                        disp('* v_r matrix:');
                        for i=1:obj.stem_data.stem_varset_r.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_r.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_r.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_r.nvar
                            for j=i+1:obj.stem_data.stem_varset_r.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_r(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_r(i,j),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end
                        disp(output);
                        output=[];
                    else
                        disp('* Pixel data are NOT cross-correlated.');
                        disp(' ');
                        disp('* Thera_r elements:');
                        output{1,1}='Variable';
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        for i=1:obj.stem_data.stem_varset_r.nvar
                            output{i+1,1}=obj.stem_data.stem_varset_r.Y_name{i};
                            output{i+1,2}=num2str(obj.stem_par.theta_r(i),'%05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{i+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                        output=[];
                    end
                end
                if obj.stem_par.k>0
                    disp(['* ',num2str(obj.stem_par.k),' fine-scale coregionalization components w_g']);
                    disp(' ');
                    disp(['* alpha_g elements:'])
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        output{i*2,1}=obj.stem_data.stem_varset_g.Y_name{i};
                        for k=1:obj.stem_par.k
                            output{i*2-1,k+1}=obj.stem_data.stem_varset_g.X_g_name{i}{k};
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i*2,k+1}=[num2str(obj.stem_par.alpha_g(i,k),'%+05.3f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f'),')'];
                            else
                                output{i*2,k+1}=num2str(obj.stem_par.alpha_g(i,k),'%+05.3f');
                            end
                            counter=counter+1;
                        end
                    end
                    disp(output);
                    output=[];
                    disp(['* theta_g elements:']);
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
                        output{k+1,2}=num2str(obj.stem_par.theta_g(k),'%06.2f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{k+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f');
                        else
                            output{k+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    output=[];
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
                        disp(['* v_g matrix for the ',num2str(k),postfix,' coreg. component:']);
                        for i=1:obj.stem_data.stem_varset_g.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_g.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_g.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_g.nvar
                            for j=i+1:obj.stem_data.stem_varset_g.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_g(i,j,k),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_g(i,j,k),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end                        
                        disp(output);
                        output=[];
                    end
                end
                if obj.stem_par.p>0
                    disp('* Transition matrix G:');
                    c=1;
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        for j=1:size(obj.stem_data.stem_varset_g.X_z{i},2)
                            output{1,c+1}=[obj.stem_data.stem_varset_g.Y_name{i},' - ',obj.stem_data.stem_varset_g.X_z_name{i}{j}];
                            output{c+1,1}=output{1,c+1};
                            c=c+1;
                        end
                    end
                    if not(isempty(obj.stem_data.stem_varset_r))
                        if not(isempty(obj.stem_data.stem_varset_r.X_z))
                            for i=1:obj.stem_data.stem_varset_r.nvar
                                for j=1:size(obj.stem_data.stem_varset_r.X_z{i},2)
                                    output{1,c+1}=[obj.stem_data.stem_varset_r.Y_name{i},' - ',obj.stem_data.stem_varset_r.X_z_name{i}{j}];
                                    output{c+1,1}=output{1,c+1};
                                    c=c+1;
                                end
                            end
                        end
                    end
                    if obj.stem_par.time_diagonal
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if i==j
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,i+1}=[num2str(obj.stem_par.G(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,i+1}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                    end
                                    counter=counter+1;
                                else
                                    output{i+1,j+1}='0';
                                end
                            end
                        end
                    else
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.G(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.G(i,j),'%+05.2f');
                                end
                                counter=counter+1;
                            end
                        end
                    end
                    disp(output);
                    output=[];
                    disp('* Sigma_eta matrix:');
                    c=1;
                    for i=1:obj.stem_data.stem_varset_g.nvar
                        for j=1:size(obj.stem_data.stem_varset_g.X_z{i},2)
                            output{1,c+1}=[obj.stem_data.stem_varset_g.Y_name{i},' - ',obj.stem_data.stem_varset_g.X_z_name{i}{j}];
                            output{c+1,1}=output{1,c+1};
                            c=c+1;
                        end
                    end
                    if not(isempty(obj.stem_data.stem_varset_r))
                        if not(isempty(obj.stem_data.stem_varset_r.X_z))
                            for i=1:obj.stem_data.stem_varset_r.nvar
                                for j=1:size(obj.stem_data.stem_varset_r.X_z{i},2)
                                    output{1,c+1}=[obj.stem_data.stem_varset_r.Y_name{i},' - ',obj.stem_data.stem_varset_r.X_z_name{i}{j}];
                                    output{c+1,1}=output{1,c+1};
                                    c=c+1;
                                end
                            end
                        end
                    end
                    if obj.stem_par.time_diagonal
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if i==j
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,i+1}=[num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,i+1}=num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f');
                                    end
                                    counter=counter+1;
                                else
                                    output{i+1,j+1}='0';
                                end
                            end
                        end
                    else
                        for j=1:obj.stem_par.p
                            for i=j:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f');
                                end
                                output{j+1,i+1}=output{i+1,j+1};
                                counter=counter+1;
                            end
                        end
                    end
                    disp(output);
                  
                end
            else
                disp('The model has not been estimated yet. Use the method print of the class stem_data to print data information.');
            end            
        end
        
        function set_system_size(obj,dim)
            %DESCRIPTION: evaluate the minimum N after which it is faster to evaluate only the diagonal of a matrix product instead of the full matrix
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the system_size property is updated
            if nargin<2
                dim=100;
                t1=0;
                t2=100;
                while t1<t2
                    m1=randn(dim);
                    m2=randn(dim);
                    tic
                    m=diag(m1*m2);
                    t1=toc;
                    tic
                    for i=1:size(m1)
                        h=m1(i,:)*m2(:,i);
                    end
                    t2=toc;
                    dim=dim+50;
                end
                obj.system_size=dim;
            else
                obj.system_size=dim;
            end
        end
        
        function set_logL(obj)
            %DESCRIPTION: evaluate the observed data log-likelihood
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the logL property of the object stem_EM_result is updated
            disp('Log-Likelihood computation...');
            st_kalman=stem_kalman(obj);
            [st_kalmansmoother_result,~,~,~,~,~,~,~,~] = st_kalman.smoother(1);
            obj.stem_EM_result.logL=st_kalmansmoother_result.logL;
            disp('Log-Likelihood computation ended.');
        end
        
        function set_varcov(obj)
            %DESCRIPTION: evaluate the variance-covariance matrix of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the varcov property of the object stem_EM_result is updated    
            
            %parameter order: beta,sigma_eps,alpha_rg,theta_r,v_r,alpha_g,theta_g,v_g,G,sigma_eta
            
            if obj.estimated==0
                error('The model has not been estimated yet');
            end
            if strcmp(obj.stem_par.correlation_type,'exponential')==0
                error('The Hessian matrix can be evaluated only in the case of exponential spatial correlation function');
            end

            N=obj.N;
            Ng=obj.Ng;
            T=obj.T;
            dim=obj.dim;            
            
            data=obj.stem_data;
            par=obj.stem_par;
            q=par.q;
            p=par.p;
            k=par.k;

            %evaluate the number of parameters
            n_beta=0;
            n_rg_alpha=0;
            n_rg_theta=0;
            n_rg_v=0;
            n_g_alpha=0;
            n_g_theta=0;
            n_g_v=0;
            n_time_G=0;
            n_time_s2e=0;
                       
            if not(isempty(data.X_rg))
                n_rg_alpha=2*q;
                if par.pixel_correlated
                    n_rg_v=q*(q-1)/2; %v_rg extra-diagonal and univoc
                    n_rg_theta=1;
                else
                    %nothing to do for V as the V matrix is the identity
                    n_rg_theta=q;
                end
            end
            
            if not(isempty(data.X_beta))
                n_beta=length(par.beta);
            end
            
            if p>0
                if par.time_diagonal
                    n_time_G=p;
                    n_time_s2e=p;
                else
                    n_time_G=p^2;
                    n_time_s2e=(p+1)*p/2;
                end
                n_time=n_time_G+n_time_s2e;
            else
                n_time=0;
            end
            
            if not(isempty(data.X_g))
                n_g_alpha=q*k; %alpha_g
                n_g_v=q*(q-1)/2*k;
                n_g_theta=k;
            end
            n_eps=size(par.sigma_eps,1);
            
            n_psi=n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+n_g_v+n_time_G+n_time_s2e;
            
            if p>0
                %kalman filter
                st_kalman=stem_kalman(obj);
                compute_logL=0;
                enable_varcov_computation=1;
                [st_kalmanfilter_result,sigma_eps,sigma_W_r,sigma_W_g,~,aj_rg,aj_g,M,sigma_geo] = st_kalman.filter(compute_logL,enable_varcov_computation);
            else
                [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,~,aj_rg,aj_g,M] = obj.get_sigma();
                st_kalmanfilter_result=stem_kalmanfilter_result([],[],[],[],[],[],[]);
            end            
            J=st_kalmanfilter_result.J(:,:,2:end); %J for t=0 is deleted
            st_kalmanfilter_result.J=st_kalmanfilter_result.J(:,:,2:end);
            st_kalmanfilter_result.zk_f=st_kalmanfilter_result.zk_f(:,2:end);
            st_kalmanfilter_result.Pk_f=st_kalmanfilter_result.Pk_f(:,:,2:end);            
           
            disp('Derivatives allocation...');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation time variant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d_P=zeros(p,p,n_psi); %the first n_beta derivatives are zero
            d_P_lag=zeros(p,p,n_psi);
            
            %d_J=zeros(p,N,n_psi); %the first n_beta derivatives are zero
            %d_J_lag=zeros(p,N,n_psi);
            for i=1:n_psi
                d_J_Lt{i}=sparse(p,N);
            end
            
            d_Z=zeros(p,n_psi); %the first n_beta derivatives are zero
            d_Z_lag=zeros(p,n_psi);
            
            d_e=zeros(N,n_psi);
            d_e_lag=zeros(N,n_psi);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            d_beta=zeros(n_beta,1,n_psi);
            d_G=zeros(p,p,n_psi);
            d_s2e=zeros(p,p,n_psi);
            
            disp('Derivatives initialization...');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  initialization of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blocks=[0 cumsum(dim)];
            %beta
            if n_beta>0
                for i=1:n_beta
                    temp=zeros(n_beta,1);
                    temp(i)=1;
                    d_beta(:,1,i)=temp;
                end
            end
            %G
            if par.time_diagonal
                for i=1:p
                    temp=zeros(p);
                    temp(i,i)=1;
                    d_G(:,:,n_psi-n_time+i)=temp;
                end
            else
                for i=1:p^2
                    temp=zeros(p^2,1);
                    temp(i)=1;
                    d_G(:,:,n_psi-n_time+i)=reshape(temp,p,p);
                end
            end
            %sigma2_eta
            if par.time_diagonal
                for i=1:p
                    temp=zeros(p);
                    temp(i,i)=1;
                    d_s2e(:,:,n_psi-n_time_s2e+i)=temp;
                end
            else
                j=1;
                for i=1:(p*(p+1))/2
                    temp=zeros(p);
                    l=1;
                    for h=1:p
                        for k=h:p
                            if i==l
                                temp(h,k)=1;
                                temp(k,h)=1;
                            end
                            l=l+1;
                        end
                    end
                    d_s2e(:,:,n_psi-n_time_s2e+i)=temp;
                    j=j+1;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %      compute d_Sgeo_prel     %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %beta
            for i=1:n_beta
                d_Sgeo_prel{i}=sparse(N,N);
            end
            %sigma_eps
            for i=n_beta+1:n_beta+n_eps
                Id=blocks(i-n_beta)+1:blocks(i-n_beta+1);
                d_Sgeo_prel{i}=sparse(Id,Id,ones(length(Id),1),N,N);
            end
            %alpha_rg
            if n_rg_alpha>0
                result=stem_misc.M_apply(sigma_W_r,M,'b');
                for j=1:n_rg_alpha
                    Id=[];
                    Jd=[];
                    elements=[];
                    for i=1:n_rg_alpha
                        temp=result(blocks(i)+1:blocks(i+1),blocks(j)+1:blocks(j+1));
                        if i==j
                            temp2=2*par.alpha_rg(j)*temp;
                            L=find(temp2);
                            [idx_I,idx_J]=ind2sub(size(temp2),L);
                            Id=[Id; idx_I+blocks(i)];
                            Jd=[Jd; idx_J+blocks(j)];
                            elements=[elements; temp2(L)];
                        else
                            temp2=par.alpha_rg(i)*temp;
                            L=find(temp2);
                            [idx_I,idx_J]=ind2sub(size(temp2),L);
                            Id=[Id; idx_I+blocks(i)];
                            Jd=[Jd; idx_J+blocks(j)];
                            elements=[elements; temp2(L)];
                            Id=[Id; idx_J+blocks(j)];
                            Jd=[Jd; idx_I+blocks(i)];
                            elements=[elements; temp2(L)];
                        end
                    end
                    zero_density=(1-length(Id)/(N^2))*100;
                    if obj.tapering||zero_density>60
                        d_Sgeo_prel{n_beta+n_eps+j}=sparse(Id,Jd,elements,N,N);
                    else
                        d_Sgeo_prel{n_beta+n_eps+j}=zeros(N);
                        temp=sub2ind([N,N],Id,Jd);
                        d_Sgeo_prel{n_beta+n_eps+j}(temp)=elements;
                    end
                end
            end
            %theta_rg
            if n_rg_theta>0
                d=stem_misc.M_apply(obj.stem_data.DistMat_r,M,'b');
                if not(par.pixel_correlated)&&(q>1)
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:2 
                            idx=blocks(j+(i-1)*2)+1:blocks(j+(i-1)*2+1);
                            result1=d(idx,idx);
                            result2=result(idx,idx);
                            result2=stem_misc.D_apply(result2,aj_rg(idx),'b');
                            result3=result1/(par.theta_r(j)^2).*result2;
                            L=find(result3);
                            [idx_I,idx_J]=ind2sub(size(result3),L);
                            Id=[Id; idx_I+blocks(j+(i-1)*2)];
                            Jd=[Jd; idx_J+blocks(j+(i-1)*2)];
                            elements=[elements; result3(L)];
                        end
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>60
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+j}=sparse(Id,Jd,elements,N,N);
                        else
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+j}=zeros(N);
                            temp=sub2ind([N,N],Id,Jd);
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+j}(temp)=elements;
                        end
                    end
                else
                    d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+1}=stem_misc.D_apply(result,aj_rg,'b').*d/(par.theta_r^2);
                end
            end
            %v_rg
            if n_rg_v>0
                if par.pixel_correlated
                    error('The variance-covariance matrix for cross-correlated pixel data is not available yet');
                    z=1;
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            for h=1:2
                                idx=blocks(j+(h-1)*2)+1:blocks(j+(h-1)*2);
                                result1=result(idx,idx);
                                result2=stem_misc.D_apply(result1,aj_rg(idx),'b')/par.v_r(i,j); 
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=[Id; idx_I+blocks(j+(h-1)*2)];
                                Jd=[Jd; idx_J+blocks(j+(h-1)*2)];
                                elements=[elements; result2(L)];
                            end
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+z}=sparse(Id,Jd,elements,N,N);
                            else
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+z}=zeros(N);
                                temp=sub2ind([N,N],Id,Jd);
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                else
                    %nothing as the matrix V_r is not estimated
                end
            end
            
            %alpha_g
            if n_g_alpha>0
                z=1;
                for k=1:par.k
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:q
                            temp=sigma_W_g{k}(blocks(i)+1:blocks(i+1),blocks(j)+1:blocks(j+1));
                            if i==j
                                temp2=2*par.alpha_g(i,k)*temp;
                                L=find(temp2);
                                [idx_I,idx_J]=ind2sub(size(temp2),L);
                                Id=[Id; idx_I+blocks(i)];
                                Jd=[Jd; idx_J+blocks(j)];
                                elements=[elements; temp2(L)];
                            else
                                temp2=par.alpha_g(i,k)*temp;
                                L=find(temp2);
                                [idx_I,idx_J]=ind2sub(size(temp2),L);
                                Id=[Id; idx_I+blocks(i)];
                                Jd=[Jd; idx_J+blocks(j)];
                                elements=[elements; temp2(L)];
                                Id=[Id; idx_J+blocks(j)];
                                Jd=[Jd; idx_I+blocks(i)];
                                elements=[elements; temp2(L)];
                            end
                        end
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>60
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=sparse(Id,Jd,elements,Ng,Ng); %note that N,N is not used
                        else
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=zeros(Ng);
                            temp=sub2ind([Ng,Ng],Id,Jd);
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}(temp)=elements;
                        end
                        z=z+1;
                    end
                end
            end
            %theta_g
            if n_g_theta>0
                for k=1:par.k
                    d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=stem_misc.D_apply(sigma_W_g{k}.*obj.stem_data.DistMat_g/(par.theta_g(k)^2),aj_g(:,k),'b');
                    if stem_misc.zero_density(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k})>60
                        d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=sparse(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k});
                    end
                end
            end
            %v_g
            if n_g_v>0
                z=1;
                for k=1:par.k
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            %since the block is extra-diagonal it has
                            %two separated D_apply!
                            result1=sigma_W_g{k}(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                            result2=stem_misc.D_apply(result1,aj_g(blocks(j)+1:blocks(j+1),k),'l');
                            result2=stem_misc.D_apply(result2,aj_g(blocks(i)+1:blocks(i+1),k),'r');
                            result2=result2/par.v_g(i,j,k);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=[Id; idx_I+blocks(j)];
                            Jd=[Jd; idx_J+blocks(i)];
                            elements=[elements; result2(L)];
                            Id=[Id; idx_J+blocks(i)];
                            Jd=[Jd; idx_I+blocks(j)];
                            elements=[elements; result2(L)];
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=sparse(Id,Jd,elements,Ng,Ng);
                            else
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=zeros(Ng);
                                temp=sub2ind([Ng,Ng],Id,Jd);
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                end
            end
            
            %G and s2e
            for i=n_psi-n_time+1:n_psi
                d_Sgeo_prel{i}=sparse(N,N);
            end
                
            if not(data.X_tv) 
                %compute d_Sgeo in the time-invariant case
                for i=1:n_beta+n_eps
                    d_Sgeo{i}=d_Sgeo_prel{i};
                end
                if n_rg_alpha>0
                    for j=1:n_rg_alpha
                        d_Sgeo{n_beta+n_eps+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+j},obj.stem_data.X_rg(:,1,1),'b');
                    end
                end
                if n_rg_theta>0
                    for j=1:n_rg_theta
                        d_Sgeo{n_beta+n_eps+n_rg_alpha+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+j},obj.stem_data.X_rg(:,1,1),'b');   
                    end
                end
                if n_rg_v>0
                    for j=1:n_rg_v
                        d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+j},obj.stem_data.X_rg(:,1,1),'b');                           
                    end
                end
                if n_g_alpha>0
                    z=1;
                    for k=1:par.k
                        for j=1:q
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z},obj.stem_data.X_g(:,1,1,k),'b');
                            L=find(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z});
                            [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}),L);
                            elements=d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}(L);
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=sparse(Id,Jd,elements,N,N);
                            z=z+1;
                        end
                    end
                end
                if n_g_theta>0
                    for k=1:par.k
                        X_g_orlated=[obj.stem_data.X_g(:,1,1,k);zeros(N-size(obj.stem_data.X_g(:,1,1,k),1),1)];
                        d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k},X_g_orlated,'b');
                    end
                end
                if n_g_v>0
                    z=1;
                    for k=1:par.k
                        for j=1:q*(q-1)/2
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z},obj.stem_data.X_g(:,1,1,k),'b');
                            L=find(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z});
                            [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}),L);
                            elements=d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}(L);
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=sparse(Id,Jd,elements,N,N);
                            z=z+1;
                        end
                    end
                end
                for i=n_psi-n_time+1:n_psi
                    d_Sgeo{i}=d_Sgeo_prel{i};   
                end
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Information matrix evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp('Hessian evaluation...');
            c0 = 0;
            IM=zeros(n_psi);
            tot=n_psi*(n_psi+1)/2*T;
            counter=1;
            
            for t=1:T
                if data.X_rg_tv
                    tRG=t;
                else
                    tRG=1;
                end
                if data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_beta_tv
                    tbeta=t;
                else
                    tbeta=1;
                end      
                if data.X_g_tv
                    tG=t;
                else
                    tG=1;
                end                  
                
                if data.X_tv
                    %compute sigma_geo in the time-variant case
                    if not(isempty(data.X_rg))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),data.X_rg(:,1,tRG),'b'),aj_rg,'b');
                    end
                    if not(isempty(data.X_g))
                        if isempty(data.X_rg)
                            if (obj.tapering)&&(p==0)
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(data.X_g,4)
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'b'),aj_g(:,k),'b');
                        end
                    end
                    if isempty(data.X_g)&&isempty(data.X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    %compute d_Sgeo in the time-variant case
                    for i=1:n_beta+n_eps
                        d_Sgeo{i}=d_Sgeo_prel{i};
                    end
                    if n_rg_alpha>0
                        for j=1:n_rg_alpha
                            d_Sgeo{n_beta+n_eps+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+j},obj.stem_data.X_rg(:,1,tRG),'b');
                        end
                    end
                    if n_rg_theta>0
                        for j=1:n_rg_theta
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+j},obj.stem_data.X_rg(:,1,tRG),'b');
                        end
                    end
                    if n_rg_v>0
                        for j=1:n_rg_v
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+j}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+j},obj.stem_data.X_rg(:,1,tRG),'b');
                        end
                    end
                    if n_g_alpha>0
                        z=1;
                        for k=1:par.k
                            for j=1:q
                                d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z},obj.stem_data.X_g(:,1,tG,k),'b');
                                L=find(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z});
                                [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}),L);
                                elements=d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}(L);
                                d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=sparse(Id,Jd,elements,N,N);
                                z=z+1;
                            end
                        end
                    end
                    if n_g_theta>0
                        for k=1:par.k
                            X_g_orlated=[obj.stem_data.X_g(:,1,tG,k);zeros(N-size(obj.stem_data.X_g(:,1,tG,k),1),1)];
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k},X_g_orlated,'b');
                        end
                    end
                    if n_g_v>0
                        z=1;
                        for k=1:par.k
                            for j=1:q*(q-1)/2
                                d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z},obj.stem_data.X_g(:,1,tG,k),'b');
                                L=find(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z});
                                [Id,Jd]=ind2sub(size(d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}),L);
                                elements=d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}(L);
                                d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=sparse(Id,Jd,elements,N,N);
                                z=z+1;
                            end
                        end
                    end
                    for i=n_psi-n_time+1:n_psi
                        d_Sgeo{i}=d_Sgeo_prel{i};
                    end
                end
                
                Lt=not(isnan(data.Y(:,t)));
                if n_time>0
                    X_z_orlated=data.X_z(:,:,tT);
                    X_z_orlated=[X_z_orlated;zeros(N-size(X_z_orlated,1),size(X_z_orlated,2))];
                    sigma_t_Lt=X_z_orlated(Lt,:)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)'+sigma_geo(Lt,Lt);
                else
                    sigma_t_Lt=sigma_geo(Lt,Lt);
                end
                
                et=zeros(N,1);
                if n_time>0
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                        et(Lt)=data.Y(Lt,t)-X_beta_orlated(Lt,:)*par.beta-X_z_orlated(Lt,:)*st_kalmanfilter_result.zk_f(:,t);
                    else
                        et(Lt)=data.Y(Lt,t)-X_z_orlated(Lt,:)*st_kalmanfilter_result.zk_f(:,t);
                    end
                else
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=[X_beta_orlated;zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2))];
                        et(Lt)=data.Y(Lt,t)-X_beta_orlated(Lt,:)*par.beta;
                    else
                        et(Lt)=data.Y(Lt,t);
                    end
                end
                
                clear sigma_eps
                clear sigma_W_r
                clear sigma_W_g
                clear result
                clear result1
                clear result2
                clear result3
                clear temp
                clear temp2
                clear Id
                clear Jd
                clear L
                clear X_g_orlated
                clear aj_g
                clear aj_rg
                clear d
                clear elements
                clear idx
                clear idx_I
                clear idx_J
                
                if t==1
                    %d_P
                    for i=n_psi-n_time_s2e+1:n_psi
                        d_P(:,:,i)=d_s2e(:,:,i);
                    end
                    
                    %d_Z is zero for each parameter at t=1
                    
                    %d_St
                    for i=1:n_psi-n_time
                        %with respect to sigma_eps
                        d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                    end
                    
                    %with respect to G are zero
                    for i=n_psi-n_time+1:n_psi-n_time_s2e
                        d_St_Lt{i}=sparse(sum(Lt),sum(Lt));
                    end
                    
                    for i=n_psi-n_time_s2e+1:n_psi
                        %with respect to sigma_eta
                        d_St_Lt{i}=X_z_orlated(Lt,:)*d_s2e(:,:,i)*X_z_orlated(Lt,:)';
                        if stem_misc.zero_density(d_St_Lt{i})>60
                            d_St_Lt{i}=sparse(d_St_Lt{i});
                        end
                    end
                    
                    %d_J
                    if n_time>0
                        for i=1:n_beta
                            d_J_Lt{i}=sparse(p,sum(Lt));
                        end
                        for i=n_beta+1:n_psi-n_time
                            d_J_Lt{i}=(-st_kalmanfilter_result.J(:,Lt,t)*d_Sgeo{i}(Lt,Lt))/sigma_t_Lt;
                        end
                        
                        for i=n_psi-n_time+1:+n_psi-n_time+n_time_G
                            %with respect to G
                            d_J_Lt{i}=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)')/sigma_t_Lt;
                        end
                        
                        for i=n_psi-n_time_s2e+1:n_psi
                            %with respect to sigma_eta
                            d_J_Lt{i}=(par.G*d_P(:,:,i)*X_z_orlated(Lt,:)')/sigma_t_Lt;
                        end
                    end
                    
                    %d_e
                    for i=1:n_beta
                        d_e(:,i)=-X_beta_orlated*d_beta(:,1,i);
                    end
                else
                    Lt1=not(isnan(data.Y(:,t-1)));
                    %d_e
                    for i=1:n_psi
                        if (i<=n_beta)
                            if n_time>0
                                d_e(:,i)=-X_beta_orlated*d_beta(:,1,i)-X_z_orlated*d_Z_lag(:,i);
                            else
                                d_e(:,i)=-X_beta_orlated*d_beta(:,1,i);
                            end
                        else
                            if n_time>0
                                d_e(:,i)=-X_z_orlated*d_Z_lag(:,i);
                            end
                        end
                    end
                    
                    %d_Z
                    if n_time>0
                        for i=1:n_beta
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_beta+1:n_psi-n_time
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time+1:n_psi-n_time+n_time_G
                            d_Z(:,i)=d_G(:,:,i)*st_kalmanfilter_result.zk_f(:,t-1)+par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time_s2e+1:n_psi
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag_Lt1{i}*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        
                        %d_P
                        for i=1:n_psi
                            d_P(:,:,i)=d_G(:,:,i)*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'-...
                                par.G*d_J_lag_Lt1{i}*X_z_orlated(Lt1,:)*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*d_P_lag(:,:,i)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*X_z_orlated(Lt1,:))*st_kalmanfilter_result.Pk_f(:,:,t-1)*d_G(:,:,i)'+d_s2e(:,:,i);
                        end
                    end
                    
                    %d_St
                    for i=1:n_psi
                        if n_time>0
                            d_St_Lt{i}=X_z_orlated(Lt,:)*d_P(:,:,i)*X_z_orlated(Lt,:)'+d_Sgeo{i}(Lt,Lt);
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                        end
                    end
                    
                    %d_J
                    if n_time>0
                        for i=1:n_psi
                            d_J_Lt{i}=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)'+par.G*d_P(:,:,i)*X_z_orlated(Lt,:)'-J(:,Lt,t)*d_St_Lt{i})/sigma_t_Lt;
                        end
                    end
                end
                
                if n_time>0
                    temp1=X_z_orlated(Lt,:)*st_kalmanfilter_result.Pk_f(:,:,t)*X_z_orlated(Lt,:)';
                    sigma_t_Lt=temp1+sigma_geo(Lt,Lt);
                else
                    sigma_t_Lt=sigma_geo(Lt,Lt);
                end
                
                clear temp0
                clear temp1
                clear temp2
                d_e_Lt=d_e(Lt,:);
                if issparse(sigma_t_Lt)
                    r = symamd(sigma_t_Lt);
                    c=chol(sigma_t_Lt(r,r));
                    temp0(r,:)=stem_misc.chol_solve(c,d_e_Lt(r,:));
                else
                    c=chol(sigma_t_Lt);
                    temp0=stem_misc.chol_solve(c,d_e_Lt);
                end
                
                for i=1:n_psi
                    temp1{i}=d_e_Lt(:,i)'*temp0;
                    d_St_i_Lt=d_St_Lt{i};
                    if nnz(d_St_i_Lt)>0
                        if issparse(sigma_t_Lt)
                            temp2{i}(r,:)=full(stem_misc.chol_solve(c,d_St_i_Lt(r,:)));
                        else
                            temp2{i}=full(stem_misc.chol_solve(c,d_St_i_Lt));
                        end
                    else
                        temp2{i}=spalloc(size(sigma_t_Lt,1),size(sigma_t_Lt,2),0);
                    end
                end
                for i=1:n_psi
                    blocks=0:50:size(temp2{i},1);
                    if not(blocks(end)==size(temp2{i},1))
                        blocks=[blocks size(temp2{i},1)];
                    end
                    nnz_temp2_i=nnz(temp2{i});
                    for j=i:n_psi
                        IM(i,j)=IM(i,j)+temp1{i}(j);
                        if (nnz_temp2_i>0)&&(nnz(temp2{j})>0)
                            sumtrace=0;
                            for z=1:length(blocks)-1
                                idx=blocks(z)+1:blocks(z+1);
                                sumtrace=sumtrace+trace(temp2{i}(idx,:)*temp2{j}(:,idx));
                            end
                            IM(i,j)=IM(i,j)+0.5*sumtrace+0.25*trace(temp2{i})*trace(temp2{j});
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
                d_P_lag=d_P;
                d_J_lag_Lt1=d_J_Lt;
                d_Z_lag=d_Z;
                d_e_lag=d_e;
                et_lag=et;
            end
            IM=IM+triu(IM,1)';
            obj.stem_EM_result.varcov=inv(IM);
        end      
        
        function set_initial_values(obj,stem_par)
            %DESCRIPTION: set the initial values of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the stem_par_initial property is updated                
            if not(isa(stem_par,'stem_par'))
                error('The input argument must be of class stem_par');
            end
            obj.stem_par_initial=stem_par;
        end
        
        %Export functions. Useful to avoid access to the properties of nested objects
        function N = N(obj)
            N=obj.stem_data.N();
        end
        
        function Nr = Nr(obj)
            Nr=obj.stem_data.Nr();
        end
        
        function Ng = Ng(obj)
            Ng=obj.stem_data.Ng();
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
        
        %Initial values estimation functions (only beta at the moment)
        function [beta0] = get_beta0(obj)
            if obj.stem_par.n_beta>0
                N = size(obj.stem_data.X_beta,1);
                y = obj.stem_data.Y(1:N,:);
                y=y(:);
                T = obj.T;
                x = zeros(N*T,size(obj.stem_data.X_beta,2));
                
                for t=1:T
                    if size(obj.stem_data.X_beta,3)==T
                        tT=t;
                    else
                        tT=1;
                    end
                    x((t-1)*N+1:t*N,:)=obj.stem_data.X_beta(:,:,tT);
                end
                
                L=not(isnan(y));
                beta0 = (x(L,:)'*x(L,:))\x(L,:)'*y(L);
            else
                disp('WARNING: the model does not include data to estimate beta');
                beta0=[];
            end
        end
      
        %Class set functions 
        function set.stem_data(obj,stem_data)
            if strcmp(class(stem_data),'stem_data')
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
            if not(length(stem_par.beta)==size(obj.stem_data.X_beta,2))
                error(['The length of beta in stem_par must be equal to ',num2str(size(obj.stem_data.X_beta,2))]);
            end
        end        
       
        function set.stem_par_initial(obj,stem_par_initial)
            if not(isa(stem_par_initial,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            else
                if obj.stem_par.time_diagonal
                    stem_par_initial.G=diag(diag(stem_par_initial.G));
                    stem_par_initial.sigma_eta=diag(diag(stem_par_initial.sigma_eta));
                end
                obj.stem_par_initial=stem_par_initial;
            end
        end
    end

end

%         function stem_par = get_initial_value_estimation(obj)
%             % beta OLS
%             if isempty(obj.stem_data.X)==0
%                 [obj.stem_par.beta res_ols] = obj.beta_ols(obj.stem_data.Y,obj.stem_data.X);
%             else
%                 res_ols=obj.stem_data.Y;
%             end
%             
%             % G and sigma_eta (only diagonal)
% %             if obj.stem_par.p>0
% %                 [G sigma_eta obj.stem_par.sigma_eps] = obj.method_of_moment(res_ols,obj.stem_data.dim);
% %                 if size(G,1)==obj.stem_par.p
% %                     obj.stem_par.G=G;
% %                     obj.stem_par.sigma_eta=sigma_eta;
% %                 else
% %                     obj.stem_par.G=mean(diag(G));
% %                     obj.stem_par.sigma_eta=mean(diag(sigma_eta));
% %                 end
% %                 
% %                 stem_kalman_obj=stem_kalman(obj);
% %                 zk_s = stem_kalman_obj.smoother();
% %                 
% %                 res_kalman = get_temporal_residual(res_ols,zk_s(:,2:end),obj.stem_par.K);
% %             else
% %                 res_kalman=res_ols;
% %             end
%  
%             % alpha and theta
% %             n_steps=30;
% %             [v theta obj.stem_par.sigma_eps] = obj.variogram(res_kalman,n_steps);
% %             obj.stem_par.alpha=sqrt(diag(v))';
% %             obj.stem_par.theta=diag(theta)';
% %             if obj.stem_par.c>0
% %                 obj.stem_par.delta=obj.stem_par.alpha/2;
% %                 obj.stem_par.theta_coreg=theta(2,1);
% %                 obj.stem_par.v=eye(obj.stem_par.q);
% %                 obj.stem_par.v(1,2)=v(2,1)/sqrt(v(1,1)*v(2,2));
% %                 obj.stem_par.v(2,1)=obj.stem_par.v(1,2);
% %             end
%             
%             %INITIAL VALUES OVERRIDING FOR TESTING PURPOSES!
%             
%             stem_par=obj.stem_par;
%             
% 
%             if obj.stem_par.p>0
%                 stem_par.G=diag(repmat(0.8,obj.stem_par.p,1));
%                 stem_par.sigma_eta=diag(repmat(0.2,obj.stem_par.p,1));
%             end
%             
%             if obj.stem_par.stem_par_constraint.sigma_eps_diag
%                 stem_par.sigma_eps=repmat(0.2,1,obj.stem_par.q);
%             else
%                 stem_par.sigma_eps=0.2*ones(obj.N,1);
%             end
%             stem_par.sigma_eps=diag(stem_par.sigma_eps);
%             if obj.stem_par.stem_par_constraint.no_direct_components==0
%                 stem_par.alpha=repmat(0.5,1,obj.stem_par.q);
%             end
%             
%             if strcmp(stem_par.correlation_type,'exponential')
%                 stem_par.theta=repmat(100,1,obj.stem_par.q);
%             end
%             if strcmp(stem_par.correlation_type,'matern')
%                 stem_par.theta=[100,0.5];
%             end
%             if stem_par.c>0
%                 stem_par.delta=repmat(0.3,1,obj.stem_par.q);
%                 if strcmp(stem_par.correlation_type,'exponential')
%                     stem_par.theta_coreg=[100]; %era 100 nella stima usata per kriging
%                 end
%                 if strcmp(stem_par.correlation_type,'matern')
%                     stem_par.theta_coreg=[30 0.5];
%                 end
%                 %v_temp=obj.stem_data.temporal_correlation('colocated');
%                 v_temp=eye(obj.stem_par.q);
%                 if min(eig(v_temp))<0
%                     [V,D]=eig(v_temp);
%                     for i=1:size(V,1)
%                         if D(i,i)<0
%                             D(i,i)=0.001;
%                         end
%                     end
%                     v_temp2=V*D*V';
%                     for j=1:size(V,1)
%                         for i=1:size(V,1)
%                             if i~=j
%                                 v_temp2(i,j)=v_temp2(i,j)/sqrt(v_temp2(i,i)*v_temp2(j,j));
%                             end
%                         end
%                     end
%                     for i=1:size(V,1)
%                         v_temp2(i,i)=1;
%                     end
%                     stem_par.v=repmat(v_temp2,[1,1,stem_par.c]);
%                     disp('Correlation matrix adjusted to be sdp');
%                     stem_par.v-v_temp
%                 else
%                     stem_par.v=repmat(v_temp,[1,1,stem_par.c]);
%                 end
%             end
%             
%             stem_par.print();
%         end
        
%         function [v,theta,sigma_eps] = variogram(obj,Y,n_steps)
%             if nargin<2
%                 Y=obj.stem_data.Y;
%             end
%             if nargin<3
%                 n_steps=50;
%             end
%             options = optimset('Display','off');
%             
%             T=obj.stem_data.T;
%             dim=obj.stem_data.dim;
%             
%             v=zeros(length(dim));
%             theta=zeros(length(dim));
%             sigma_eps=zeros(length(dim));
%             
%             blocks=[0 cumsum(dim)];
%             for k=1:length(dim)
%                 for h=k:length(dim)
%                     y_sub1=Y(blocks(h)+1:blocks(h+1),:); 
%                     y_sub2=Y(blocks(k)+1:blocks(k+1),:);
%                     DistMat_sub=obj.stem_data.DistMat(blocks(h)+1:blocks(h+1),blocks(k)+1:blocks(k+1));
%                     max_distance=max(obj.stem_data.DistMat(:));
%                     step=max_distance/n_steps;
%                     groups=0:step:max_distance;
%                     group=zeros(length(groups)-1,length(dim));
%                     count=zeros(length(groups)-1,length(dim));
%                     
%                     index_m=1;
%                     for t=1:T
%                         if (sum(isnan(y_sub1(:,t)))/length(y_sub1(:,t))<0.8)&&(sum(isnan(y_sub2(:,t)))/length(y_sub2(:,t))<0.8)
%                             index=1;
%                             tot=dim(k)*dim(h);
%                             vx=zeros(tot,1);
%                             vy=zeros(tot,1);
%                             for j=1:dim(k)
%                                 for i=1:dim(h)
%                                     vx(index)=DistMat_sub(i,j);
%                                     vy(index)=0.5*(y_sub1(i,t)-y_sub2(j,t))^2;
%                                     index=index+1;
%                                 end
%                             end
%                             for i=2:length(groups)
%                                 data=vy(vx>=groups(i-1)&vx<groups(i));
%                                 group(i-1,k)=nanmean(data);
%                                 count(i-1,k)=length(data(isnotnan(data)));
%                             end
%                             xv=(groups(2:8)+(step/2))';
%                             yv=group(2:8,k);
%                             xv=xv(not(isnan(yv)));
%                             yv=yv(not(isnan(yv)));
%                             b=regress(yv,[ones(length(xv),1) xv]);
%                             v0(index_m)=b(1);
%                             
%                             L=not(isnan(group(:,k)));
%                             L(1)=0;
%                             g=groups(1:end-1)+(step/2);
%                             x1=g(L)';
%                             y1=group(L,k);
%                             weights=count(L,k);
%                             x0=[1 20];
%                             
%                             if 1
%                                 f = @(x,xdata) x(1)*(1-exp(-xdata/x(2)))+v0(index_m);'x';'xdata';
%                                 c = lsqcurvefit(f,x0,x1,y1,[],[],options);
%                                 theta0(index_m)=c(2);
%                                 vv(index_m)=c(1);
%                                 index_m=index_m+1;
%                             else
%                                 o = fitoptions('Method','NonlinearLeastSquares');
%                                 o.Weights=weights;
%                                 o.StartPoint=[1 20];
%                                 f = fittype('a*(1-exp(-x/b))');
%                                 try
%                                     fit1 = fit(x1,y1-v0(index_m),f,o);
%                                     c=coeffvalues(fit1);
%                                     theta0(index_m)=c(2);
%                                     vv(index_m)=c(1);
%                                     index_m=index_m+1;
%                                 catch
%                                     disp('Error');
%                                 end
%                             end
%                             
%                             
%                             if 0
%                                 xx=0:0.1:max(groups);
%                                 yy=c(1)*(1-exp(-xx./c(2)))+v0(index_m);
%                                 plot(groups(1:end-1)+(step/2),group(:,k),'o');
%                                 hold on
%                                 plot(xx,yy,'-r');
%                             end
%                         end
%                     end
%                     filter=(theta0<max_distance/0.75)&(v0>0);
%                     theta0=theta0(filter);
%                     vv=vv(filter); %filtro anche v_tot perch� stimato congiuntamente a theta
%                     v0=v0(filter);
%                     
%                     sigma_eps(h,k)=mean(v0);
%                     v(h,k)=mean(vv);
%                     theta(h,k)=mean(theta0);   
%                 end
%             end
%         end

%         function data_reset(obj)
%             obj.stem_data.reset();
%         end
        
%         function plot_par(obj,iterations)
%             if isempty(obj.stem_par_all)==0
%                 n_sub=2*obj.stem_par.q+length(obj.stem_par.beta)+2*obj.stem_par.p^2;
%                 if strcmp(obj.stem_par.correlation_type,'exponential')
%                     n_sub=n_sub+obj.stem_par.q;
%                 end
%                 if strcmp(obj.stem_par.correlation_type,'matern')
%                     n_sub=n_sub+obj.stem_par.q*2;
%                 end
%                 if obj.stem_data.simulated
%                     n_sub=n_sub+obj.stem_par.q; %for the error term
%                 end
%                 if obj.stem_par.c>0
%                     n_sub=n_sub+obj.stem_par.q+(obj.stem_par.q*(obj.stem_par.q-1))/2;
%                     if strcmp(obj.stem_par.correlation_type,'exponential')
%                         n_sub=n_sub+obj.stem_par.c;
%                     end
%                     if strcmp(obj.stem_par.correlation_type,'matern')
%                         n_sub=n_sub+2*obj.stem_par.c;
%                     end
%                     if obj.stem_data.simulated
%                         n_sub=n_sub+obj.stem_par.c;%for the error term
%                     end
%                 end
%                 l=n_sub^0.5;
%                 if round(l)^2==n_sub
%                     rows=l;
%                     cols=l;
%                 else
%                     rows=ceil(l);
%                     cols=round(l);
%                 end
%                 
%                 counter=1;
%                 for i=1:length(obj.stem_par.beta)
%                     subplot(rows,cols,counter);
%                     plot(obj.stem_par_all(counter,:));
%                     title(['beta ',num2str(counter)]);
%                     xlim([1 iterations]);
%                     counter=counter+1;
%                 end
%                 for i=1:obj.stem_par.q
%                     subplot(rows,cols,counter);
%                     plot(obj.stem_par_all(counter,:));
%                     title(['sigma eps ',obj.stem_data.variable_name{i}]);
%                     xlim([1 iterations]);
%                     counter=counter+1;
%                 end
%                 for i=1:obj.stem_par.q
%                     subplot(rows,cols,counter);
%                     plot(obj.stem_par_all(counter,:));
%                     title(['alpha ',obj.stem_data.variable_name{i}]);
%                     xlim([1 iterations]);
%                     counter=counter+1;
%                 end
%                 for i=1:obj.stem_par.q
%                     if strcmp(obj.stem_par.correlation_type,'exponential')
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['theta ',obj.stem_data.variable_name{i}]);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                     if strcmp(obj.stem_par.correlation_type,'matern')
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['matern alpha ',obj.stem_data.variable_name{i}]);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                         
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['matern nu ',obj.stem_data.variable_name{i}]);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                 end
%                 if obj.stem_par.c>0
%                     for i=1:obj.stem_par.q
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['delta ',obj.stem_data.variable_name{i}]);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                     for i=1:obj.stem_par.c
%                         if strcmp(obj.stem_par.correlation_type,'exponential')
%                             subplot(rows,cols,counter);
%                             plot(obj.stem_par_all(counter,:));
%                             title(['theta coreg ',num2str(i)]);
%                             xlim([1 iterations]);
%                             counter=counter+1;
%                         end
%                         if strcmp(obj.stem_par.correlation_type,'matern')
%                             subplot(rows,cols,counter);
%                             plot(obj.stem_par_all(counter,:));
%                             title(['matern coreg alpha ',num2str(i)]);
%                             xlim([1 iterations]);
%                             counter=counter+1;
%                             subplot(rows,cols,counter);
%                             plot(obj.stem_par_all(counter,:));
%                             title(['matern coreg nu ',num2str(i)]);
%                             xlim([1 iterations]);
%                             counter=counter+1;
%                         end
%                     end
%                     for i=1:(obj.stem_par.q*(obj.stem_par.q-1))/2
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['v ',num2str(i)]);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                 end
%                 if obj.stem_par.p>0
%                     for i=1:obj.stem_par.p^2
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['G ']);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                     for i=1:obj.stem_par.p^2
%                         subplot(rows,cols,counter);
%                         plot(obj.stem_par_all(counter,:));
%                         title(['sigma eta ']);
%                         xlim([1 iterations]);
%                         counter=counter+1;
%                     end
%                 end
%                 drawnow
%             end
%         end  %dovrebbe diventare un metodo di stem_par?

%     methods (Static)
%         
%         function [G sigma_eta sigma_eps]=method_of_moment(Y,dim)
%             n_var=length(dim);
%             
%             row_start=1;
%             for j=1:n_var
%                 index=1;
%                 for i=row_start:row_start+dim(j)-1
%                     if sum(isnan(Y(i,:)))/length(Y(i,:))<0.9 %se i missing sono meno del 50% allora stimo g
%                         L=not(isnan(Y(i,:)));
%                         data=Y(i,L);
%                         if length(data)>=1
%                             c=xcorr(data,2,'unbiased');
%                             g0(index)=c(end)/c(end-1);
%                             temp=xcov(data,1,'unbiased');
%                             gamma0(index)=temp(end-1);
%                             gamma1(index)=temp(end);
%                             index=index+1;
%                         end
%                     end
%                 end
%                 filter=abs(g0)<1;
%                 g0=g0(filter);
%                 gamma0=gamma0(filter);
%                 gamma1=gamma1(filter);
%                 for i=1:length(g0)
%                     s2eta0(i)=gamma1(i)*(1-g0(i)^2)/g0(i);
%                     s2eps0(i)=gamma0(i)-(s2eta0(i)/(1-g0(i)^2));
%                 end
%                 s2eps0=s2eps0(s2eps0>0);
%                 s2eps0=median(s2eps0);
%                 s2eta0=s2eta0(s2eta0>0);
%                 s2eta0=median(s2eta0);
%                 g0=median(g0);
%                 
%                 G(j,j)=g0;
%                 sigma_eta(j,j)=s2eta0;
%                 sigma_eps(j,j)=s2eps0;
%                 row_start=row_start+dim(j)-1;
%             end
%         end
%     end

