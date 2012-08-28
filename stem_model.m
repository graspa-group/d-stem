%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_model < handle
    
    %CONSTANTS
    %N  = n1_g+...+nq_g+n1_r+...+nq_r - total number of observation sites
    %Ng = n1_g+...+nq_g - total number of ground level sites
    %Nr = n1_r+...+nq_r - total number of remote sensing sites
    %M = 2 if both ground level and remote sense data are considered. M = 1 if only ground level data are considered.
    
    properties
        stem_data=[];           %[stem_data object] (1x1) object containing all the data used to estimated the model
        stem_par=[];            %[stem_par object]  (1x1) updated at each iteration of the EM algorithm
        stem_par_initial=[];    %[stem_par object]  (1x1) parameter starting values
        stem_par_sim=[];        %[stem_par object]  (1x1) parameter values used to simulate data
        note=[];
                estimated=0;            %[boolean] (1x1) 0: the model is not estimated; 
    end
    
    properties (SetAccess = private)
        stem_EM_result=[];      %[stem_EM_result object] (1x1) object containing all the results of the EM estimation

                                                %1: the model has been estimated.
        cross_validation=0;     %[boolean] (1x1) 0: the model has been estimated considering all the data; 
                                                %1: the model has bee estimated excluding the cross-validation data.
        system_size=100;        %[integer] (1x1) if N > system_size than only the diagonal is computed in matrix multiply operations
        tapering=[];
        tapering_r=[];
        tapering_g=[];
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
        end
        
        function [aj_rg,aj_g] = get_aj(obj)
            %DESCRIPTION: provides the vector aj_rg and aj_g useful in the EM estimation
            %
            %INPUT
            %obj   - [stem_model object] (1x1)
            %
            %OUTPUT
            %aj_rg - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_rg
            %aj_g  - [double]            (Nx1) is the diagonal of the NxN diagonal matrix a*J_g. 
            
            %NOTE
            %The elements of aj_g from Ng+1 to N are all zeros. This allows the
            %use of stem_misc.D_apply both for the remote sensing data and the ground
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
            %DESCRIPTION: provides the vectors aj_rg_r and j_r useful in the EM estimation
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
            %DESCRIPTION: provides the vectors aj_g_rs and j_r useful in the EM estimation
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
            %DESCRIPTION: provides the variance-covariance matrices and some 
            %vectors that are needed by the EM algorithm
            %
            %INPUT
            %
            %obj         - [stem_model object] (1x1)
            %<sigma_W_r> - [double]            (NrXNr) (default: []) variance-covariance matrix of W_r. It is provided as input argument during kriging when sigma_W_r does not change across blocks  
            %
            %OUTPUT
            %
            %sigma_eps   - [double]            (NxN) variance-covariance matrix of epsilon
            %sigma_W_r   - [double]            (NrxNr) variance-covariance matrix of W_r
            %sigma_W_g   - [double]            (NgxN1xq) variance-covariance matrices of the q W_g_i
            %sigma_geo   - [double]            (NxN) variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %sigma_Z     - [double]            (pxp) variance-covariance of Z
            %aj_rg       - [double]            (Nx1) see the details of the function get_aj;
            %aj_g        - [double]            (Nx1) see the details of the function get_aj;
            %M           - [integer >0]        (Ngx1) see the details of the function update_M of the class stem_data
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
            
            %CONTROLLARE LA CREAZIONE DI SIGMA_W_R NEL CASO REMOTE_CORRELATED = FALSE!!!
            
            %sigma_W_r
            if not(isempty(obj.stem_data.stem_varset_r))
                if not(isempty(obj.stem_data.stem_varset_r.X_rg))
                    if (nargin==1)||isempty(sigma_W_r)
                        if not(obj.tapering_r)
                            sigma_W_r=zeros(obj.stem_data.stem_varset_r.N);
                        end
                        blocks=[0 cumsum(obj.stem_data.stem_varset_r.dim)];
                        if obj.stem_par.remote_correlated
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
                                        size=length(corr_result.I);
                                        I(idx+1:idx+size)=corr_result.I+blocks(i);
                                        J(idx+1:idx+size)=corr_result.J+blocks(j);
                                        elements(idx+1:idx+size)=corr_result.correlation;
                                        idx=idx+size;
                                        if not(i==j)
                                            I(idx+1:idx+size)=corr_result.J+blocks(j);
                                            J(idx+1:idx+size)=corr_result.I+blocks(i);
                                            elements(idx+1:idx+size)=corr_result.correlation;
                                            idx=idx+size;
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
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_r,obj.stem_data.DistMat_r(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                    weights=stem_misc.wendland(obj.stem_data.DistMat_r(idx_rc,idx_rc),obj.stem_data.stem_gridlist_r.tap);
                                    corr_result.correlation=obj.stem_par.v_r(i,i)*corr_result.correlation.*weights;
                                    size=length(corr_result.I);
                                    I(idx+1:idx+size)=corr_result.I+blocks(i);
                                    J(idx+1:idx+size)=corr_result.J+blocks(i);
                                    elements(idx+1:idx+size)=corr_result.correlation;
                                    idx=idx+size;
                                end
                            end
                            if obj.tapering_r
                                sigma_W_r=sparse(I,J,elements);
                            end
                        end
                    end
                end
            end
            
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
                               size=length(corr_result.I);
                               I(idx+1:idx+size)=corr_result.I+blocks(i);
                               J(idx+1:idx+size)=corr_result.J+blocks(j);
                               elements(idx+1:idx+size)=corr_result.correlation;
                               idx=idx+size;
                               if not(i==j)
                                   I(idx+1:idx+size)=corr_result.J+blocks(j);
                                   J(idx+1:idx+size)=corr_result.I+blocks(i);
                                   elements(idx+1:idx+size)=corr_result.correlation;
                                   idx=idx+size;
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
           
           %sigma_geo
           [aj_rg,aj_g]=obj.get_aj;
           if not(obj.stem_data.X_tv)
               %time invariant case
               if not(isempty(obj.stem_data.stem_varset_r))
                   sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),obj.stem_data.X_rg(:,1,1),'b'),aj_rg,'b');
               end
               if obj.stem_par.k>0
                   if isempty(obj.stem_data.stem_varset_r)
                       %se manca il remoto allora sigma_geo non è stata
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
            %<nan_rate>        - [double [0,1]]      (Mx1) (default: []) missing rate of the simulated data
            %<nan_pattern_par> - [double >0]         (Mx1) (default: []) (UoM: km) parameter of the exponential spatial correlation function used to define the spatial pattern of the missing data
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
            %obj                - [stem_model object] (1x1)
            %<exit_toll>        - [double >0]         (1x1) (default: 0.0001) the EM algorithm stops if the relative norm between two consecutive iterations is below exit_toll
            %<max_iterations>   - [integer >0]        (1x1) (default: 1000)  the EM algorithm stops if the number of iterations exceed max_iterations
            %<numeric_opt_type> - [string]            (1x1) (default: 'single') if 'single' then elements of the V_i matrices are numerically estimated one-by-one. If 'full' the elements are jointly estimated.
            %<pathparallel>     - [string]            (1x1) (default: []) full or relative path of the folder to use for parallel computation
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
                
                %recover the cross-validation sites from the step defined
                %in the stem_crossval object
                step=obj.stem_data.stem_crossval.idx_step;
                indices=1:size(obj.stem_data.stem_varset_g.Y{idx_var},1);
                indices(not(mod(indices,step)==0))=[];
                
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
                if not(isempty(obj.stem_data.stem_varset_g.X_time))
                    X_time={obj.stem_data.stem_varset_g.X_time{idx_var}(indices,:,:)};
                    X_time_name={obj.stem_data.stem_varset_g.X_time_name{idx_var}};
                else
                    X_time={};
                    X_time_name={};
                end          
                if not(isempty(obj.stem_data.stem_varset_g.X_g))
                    X_g={obj.stem_data.stem_varset_g.X_g{idx_var}(indices,:,:,:)};
                    X_g_name={obj.stem_data.stem_varset_g.X_g_name{idx_var}};
                else
                    X_g={};
                    X_g_name={};
                end      
                
                %set the cross_mindistance vector
                dim=obj.stem_data.dim;
                blocks=[0 cumsum(dim)];
                temp_dist=obj.stem_data.DistMat_g(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                temp_dist=temp_dist(indices,:);
                temp_dist(:,indices)=[];
                obj.stem_data.stem_crossval.mindistance=min(temp_dist');
                clear temp_dist
                
                obj.stem_data.stem_crossval.stem_varset=stem_varset(Y,Y_name,X_rg,X_rg_name,X_beta,X_beta_name,X_time,X_time_name,X_g,X_g_name);
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
            %set the current parameter value with the estimated initial
            %value
            obj.stem_par=obj.stem_par_initial;
            if isempty(stem_EM_options.pathparallel)
                obj.stem_EM_result=st_EM.estimate();
            else
                obj.stem_EM_result=st_EM.estimate_parallel(stem_EM_options.pathparallel);
            end
            obj.estimated=1;
            if obj.cross_validation
                st_krig=stem_krig(obj);
                block_size=1100;
                back_transform=0;
                no_varcov=1;
                crossval=1;
                obj.stem_data.stem_crossval.stem_krig_result=st_krig.kriging(obj.stem_data.stem_crossval.variable_name,[],block_size,[],[],back_transform,no_varcov,crossval);
                res=obj.stem_data.stem_crossval.stem_krig_result.y_hat-obj.stem_data.stem_crossval.stem_varset.Y{1};
                obj.stem_data.stem_crossval.mse=nanvar(res(:));
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
            disp('Log-Likelihood computation...');
            st_kalman=stem_kalman(obj);
            [st_kalmansmoother_result,~,~,~,~,~,~,~,~] = st_kalman.smoother(1);
            obj.stem_EM_result.logL=st_kalmansmoother_result.logL;
            disp('Log-Likelihood computation ended.');
        end
        
        function set_varcov(obj)
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
                if par.remote_correlated
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
                n_g_v=q*(q-1)*k;
                n_g_theta=k;
            end
            n_eps=size(par.sigma_eps,1);
            
            n_psi=n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+n_g_v+n_time_G+n_time_s2e;
            
            if p>0
                %kalman filter
                st_kalman=stem_kalman(obj);
                [st_kalmanfilter_result,sigma_eps,sigma_W_r,sigma_W_g,~,aj_rg,aj_g,M,sigma_geo] = st_kalman.filter();
            else
                [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,~,aj_rg,aj_g,M] = obj.get_sigma();
                st_kalmanfilter_result=stem_kalmanfilter_result([],[],[],[],[],[]);
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
            
            d_J=zeros(p,N,n_psi); %the first n_beta derivatives are zero
            d_J_lag=zeros(p,N,n_psi);
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %  compute d_Sgeo_prel  %
            %%%%%%%%%%%%%%%%%%%%%%%%%

            %beta
            for i=1:n_beta
                if obj.tapering
                    d_Sgeo_prel{i}=sparse(N,N); 
                else
                    d_Sgeo_prel{i}=zeros(N); 
                end
            end
            %sigma_eps
            for i=n_beta+1:n_beta+n_eps
                if obj.tapering
                    Id=blocks(i-n_beta)+1:blocks(i-n_beta+1);
                    d_Sgeo_prel{i}=sparse(Id,Id,ones(length(Id),1),N,N);
                else
                    d_Sgeo_prel{i}=zeros(N);
                    d_Sgeo_prel{i}(blocks(i-n_beta)+1:blocks(i-n_beta+1),blocks(i-n_beta)+1:blocks(i-n_beta+1))=eye(dim(i-n_beta));
                end
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
                    if obj.tapering
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
                if not(par.remote_correlated)
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:2 
                            idx=blocks(j*i)+1:blocks(j*i+1);
                            result1=d(idx,idx);
                            result2=result(idx,idx);
                            result2=stem_misc.D_apply(result2,aj_rg(idx),'b');
                            result3=result1/(par.theta_r(j)^2).*result2;
                            L=find(result3);
                            [idx_I,idx_J]=ind2sub(size(result3),L);
                            Id=[Id; idx_I+blocks(j*i)];
                            Jd=[Jd; idx_J+blocks(j*i)];
                            elements=[elements; result3(L)];
                        end
                        if obj.tapering
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
                if par.remote_correlated
                    z=1;
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            for h=1:2
                                idx=blocks(j*h)+1:blocks(j*h+1);
                                result1=result(idx,idx);
                                result2=stem_misc.D_apply(result1,aj_rg(idx),'b')/par.v_r(i,j);
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=[Id; idx_I+blocks(j*h)];
                                Jd=[Jd; idx_J+blocks(j*h)];
                                elements=[elements; result2(L)];
                            end
                            if obj.tapering
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
                                temp2=2*par.alpha_g(j,k)*temp;
                                L=find(temp2);
                                [idx_I,idx_J]=ind2sub(size(temp2),L);
                                Id=[Id; idx_I+blocks(i)];
                                Jd=[Jd; idx_J+blocks(j)];
                                elements=[elements; temp2(L)];
                            else
                                temp2=alpha_g(j,k)*temp;
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
                        if obj.tapering
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=sparse(Id,Jd,elements,Ng,Ng); %note that N,N is not used
                        else
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=zeros(Ng);
                            elements=reshape(elements,Ng,Ng);
                            d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+z}=elements;
                        end
                        z=z+1;
                    end
                end
            end
            %theta_g
            if n_g_theta>0
                for k=1:par.k
                    d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=sigma_W_g{k}.*obj.stem_data.DistMat_g/(par.theta_g(k)^2);
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
                            for h=1:2
                                result1=sigma_W_g{k}(blocks(j*h-1)+1:blocks((j*h+1)-1),blocks(i*h-1)+1:blocks((i*h+1)-1));
                                result2=stem_misc.D_apply(result1,aj_g(:,k))/par.v_g(i,j,k);
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=[Id; idx_I+blocks(j*h)];
                                Jd=[Jd; idx_J+blocks(j*h)];
                                elements=[elements; result2(L)];
                            end
                            if obj.tapering
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=sparse(Id,Jd,elements,Ng,Ng);
                            else
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=zeros(Ng);
                                elements=reshape(elements,length(Id),length(Jd));
                                d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=elements;
                            end
                            z=z+1;
                        end
                    end
                end
            end
            
            %G and s2e
            for i=n_psi-n_time+1:n_psi
                if obj.tapering
                    d_Sgeo_prel{i}=sparse(N,N);
                else
                    d_Sgeo_prel{i}=zeros(N);
                end
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
                        for j=1:p
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
                        d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=stem_misc.D_apply(stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k},obj.stem_data.X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                    end
                end
                if n_g_v>0
                    z=1;
                    for k=1:par.k
                        for j=1:q*(q-1)
                            d_Sgeo{n_beta+n_eps+n_rg_alpha_r+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha_r+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta},obj.stem_data.X_g(:,1,1,k),'b');
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
                if data.X_time_tv
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
                            d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k}=stem_misc.D_apply(stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+k},obj.stem_data.X_g(:,1,tG,k),'b'),aj_g(:,k),'b');
                        end
                    end
                    if n_g_v>0
                        z=1;
                        for k=1:par.k
                            for j=1:q*(q-1)
                                d_Sgeo{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z}=stem_misc.D_apply(d_Sgeo_prel{n_beta+n_eps+n_rg_alpha+n_rg_theta+n_rg_v+n_g_alpha+n_g_theta+z},obj.stem_data.X_g(:,1,tG,k),'b');
                            end
                        end
                    end
                    for i=n_psi-n_time+1:n_psi
                        d_Sgeo{i}=d_Sgeo_prel{i};
                    end
                end
                
                if n_time>0
                    sigma_t=data.X_time(:,:,tT)*st_kalmanfilter_result.Pk_f(:,:,t)*data.X_time(:,:,tT)'+sigma_geo;
                else
                    sigma_t=sigma_geo;
                end
                et=zeros(N,1);
                Lt=not(isnan(data.Y(:,t)));
                
                if n_time>0
                    if n_beta>0
                        et(Lt)=data.Y(Lt,t)-data.X_beta(Lt,:,tbeta)*par.beta-data.X_time(Lt,:,tT)*st_kalmanfilter_result.zk_f(:,t);
                    else
                        et(Lt)=data.Y(Lt,t)-data.X_time(Lt,:,tT)*st_kalmanfilter_result.zk_f(:,t);
                    end
                else
                    if n_beta>0
                        et(Lt)=data.Y(Lt,t)-data.X_beta(Lt,:,tbeta)*par.beta;
                    else
                        et(Lt)=data.Y(Lt,t);
                    end
                end
                
                if t==1
                    %d_P
                    for i=n_psi-n_time_s2e+1:n_psi
                        d_P(:,:,i)=d_s2e(:,:,i);
                    end
                    
                    %d_Z is zero for each parameter at t=1
                    
                    %d_St
                    for i=1:n_psi-n_time
                        %with respect to sigma_eps
                        d_St{i}=d_Sgeo{i};
                    end
                    
                    %with respect to G are zero
                    for i=n_psi-n_time+1:n_psi-n_time_s2e
                        if obj.tapering
                            d_St{i}=sparse(N,N);
                        else
                            d_St{i}=zeros(N);
                        end
                    end
                    
                    for i=n_psi-n_time_s2e+1:n_psi
                        %with respect to sigma_eta
                        d_St{i}=data.X_time(:,:,tT)*d_s2e(:,:,i)*data.X_time(:,:,tT)';
                    end
                    
                    %d_J
                    if n_time>0
                        for i=n_beta+1:n_psi-n_time
                            d_J(:,Lt,i)=(-st_kalmanfilter_result.J(:,Lt,t)*d_Sgeo{i}(Lt,Lt))/sigma_t(Lt,Lt);
                        end
                        
                        for i=n_psi-n_time+1:+n_psi-n_time+n_time_G
                            %with respect to G
                            d_J(:,Lt,i)=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*data.X_time(Lt,:,tT)')/sigma_t(Lt,Lt);
                        end
                        
                        
                        for i=n_psi-n_time_s2e+1:n_psi
                            %with respect to sigma_eta
                            d_J(:,Lt,i)=(par.G*d_P(:,:,i)*data.X_time(Lt,:,tT)')/sigma_t(Lt,Lt);
                        end
                    end
                    
                    %d_e
                    for i=1:n_beta
                        d_e(:,i)=-data.X_beta(:,:,tbeta)*d_beta(:,1,i);
                    end
                else
                    Lt1=not(isnan(data.Y(:,t-1)));
                    %d_e
                    for i=1:n_psi
                        if (i<=n_beta)
                            if n_time>0
                                d_e(:,i)=-data.X_beta(:,:,tbeta)*d_beta(:,1,i)-data.X_time(:,:,tT)*d_Z_lag(:,i);
                            else
                                d_e(:,i)=-data.X_beta(:,:,tbeta)*d_beta(:,1,i);
                            end
                        else
                            if n_time>0
                                d_e(:,i)=-data.X_time(:,:,tT)*d_Z_lag(:,i);
                            end
                        end
                    end
                    
                    %d_Z
                    if n_time>0
                        for i=1:n_beta
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_beta+1:n_psi-n_time
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag(:,Lt1,i)*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time+1:n_psi-n_time+n_time_G
                            d_Z(:,i)=d_G(:,:,i)*st_kalmanfilter_result.zk_f(:,t-1)+par.G*d_Z_lag(:,i)+d_J_lag(:,Lt1,i)*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        for i=n_psi-n_time_s2e+1:n_psi
                            d_Z(:,i)=par.G*d_Z_lag(:,i)+d_J_lag(:,Lt1,i)*et_lag(Lt1)+J(:,Lt1,t-1)*d_e_lag(Lt1,i);
                        end
                        
                        %d_P
                        for i=1:n_psi
                            d_P(:,:,i)=d_G(:,:,i)*(eye(p)-J(:,Lt1,t-1)*data.X_time(Lt1,:,tT))*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'-...
                                par.G*d_J_lag(:,Lt1,i)*data.X_time(Lt1,:,tT)*st_kalmanfilter_result.Pk_f(:,:,t-1)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*data.X_time(Lt1,:,tT))*d_P_lag(:,:,i)*par.G'+...
                                par.G*(eye(p)-J(:,Lt1,t-1)*data.X_time(Lt1,:,tT))*st_kalmanfilter_result.Pk_f(:,:,t-1)*d_G(:,:,i)'+d_s2e(:,:,i);
                        end
                    end
                    
                    %d_St
                    for i=1:n_psi
                        if n_time>0
                            d_St{i}(Lt,Lt)=data.X_time(Lt,:,tT)*d_P(:,:,i)*data.X_time(Lt,:,tT)'+d_Sgeo{i}(Lt,Lt);
                        else
                            d_St{i}(Lt,Lt)=d_Sgeo{i}(Lt,Lt);
                        end
                    end
                    
                    %d_J
                    if n_time>0
                        for i=1:n_psi
                            d_J(:,Lt,i)=(d_G(:,:,i)*st_kalmanfilter_result.Pk_f(:,:,t)*data.X_time(Lt,:,tT)'+par.G*d_P(:,:,i)*data.X_time(Lt,:,tT)'-J(:,Lt,t)*d_St{i}(Lt,Lt))/sigma_t(Lt,Lt);
                        end
                    end
                end
                
                if n_time>0
                    if data.X_time_tv
                        temp1=data.X_time(:,:,t)*st_kalmanfilter_result.Pk_f(:,:,t)*data.X_time(:,:,t)';
                    else
                        temp1=data.X_time*st_kalmanfilter_result.Pk_f(:,:,t)*data.X_time';
                    end
                    sigma_t=temp1+sigma_geo;
                else
                    sigma_t=sigma_geo;
                end

                S=sigma_t(Lt,Lt);
                r = symamd(S);
                chol_S=chol(S(r,r));
                
                clear temp0
                clear temp1
                clear temp2
                d_e_Lt=d_e(Lt,:);
                temp0(r,:)=stem_misc.chol_solve(chol_S,d_e_Lt(r,:));
                for i=1:n_psi
                    tic
                    temp1{i}=d_e_Lt(:,i)'*temp0;
                    d_St_i_Lt=d_St{i}(Lt,Lt);
                    temp2{i}(r,:)=stem_misc.chol_solve(chol_S,d_St_i_Lt(r,r),1);
                    toc
                end
                for i=1:n_psi
                    for j=i:n_psi
                        IM(i,j)=IM(i,j)+temp1{i}(j);
                        IM(i,j)=IM(i,j)+0.5*trace(temp2{i}*temp2{j})+0.25*trace(temp2{i})*trace(temp2{j});
                        counter=counter+1;
                        if (mod(counter,100)==0)||(counter==tot)
                            if ((mod(round(counter/tot*100),20)==0)||(round(counter/tot*100)<3)) ...
                                && c0 ~= round(counter/tot*100);
                                c0 = round(counter/tot*100);
                                disp(['Hessian evaluation: ',num2str(round(counter/tot*100)),'% completed']);
                            end
                        end
                    end
                end
                d_P_lag=d_P;
                d_J_lag=d_J;
                d_Z_lag=d_Z;
                d_e_lag=d_e;
                et_lag=et;
            end
            IM=IM+triu(IM,1)';
            obj.stem_EM_result.varcov=inv(IM);
        end      
        
        function set_initial_values(obj,stem_par)
            if not(isa(stem_par,'stem_par'))
                error('The input argument must be of class stem_par');
            end
            obj.stem_par_initial=stem_par;
        end
        
        %export functions
        function N = N(obj)
            %return the value of N from the stem_data object
            N=obj.stem_data.N();
        end
        
        function Nr = Nr(obj)
            Nr=obj.stem_data.Nr();
        end
        
        function Ng = Ng(obj)
            Ng=obj.stem_data.Ng();
        end
        
        function T = T(obj)
            %return the value of T from the stem_data object
            T=obj.stem_data.T();
        end
        
        function nvar=nvar(obj)
            %return the value of nvar from the stem_data object
            nvar=obj.stem_data.nvar();
        end
        
        function dim=dim(obj)
            %return the value of dim from the stem_data object
            dim=obj.stem_data.dim();
        end
        
        %initial values estimation functions
        
        function [beta0] = get_beta0(obj)
            N = obj.N;
            y = obj.stem_data.Y(:);
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
%                     vv=vv(filter); %filtro anche v_tot perchè stimato congiuntamente a theta
%                     v0=v0(filter);
%                     
%                     sigma_eps(h,k)=mean(v0);
%                     v(h,k)=mean(vv);
%                     theta(h,k)=mean(theta0);   
%                 end
%             end
%         end

%         function st_loo_residual = leave_one_out(obj,stem_par_alldata,EM_toll,EM_iter,varset)
%             if nargin<3
%                 EM_iter=100;
%             end
%             if nargin<4
%                 EM_toll=0.001;
%             end
%             if EM_iter<=0
%                 error('EM_iter must be >0');
%             end
%             if EM_toll<=0
%                 error('EM_toll must be >0');
%             end
%             if (nargin<2)||isempty(stem_par_alldata)
%                 cdisp('Full data model estimation started...');
%                 obj.stem_par_initial=obj.get_initial_value_estimation();
%                 obj.EM_estimate(EM_toll,EM_iter,'single');
%                 kriging_Var_W_bar_hat=obj.stem_EM_result.Var_W_bar_hat;
%                 cdisp('Full data model estimation ended');
%                 stem_par_alldata=obj.stem_par;
%             else
%                 kriging_Var_W_bar_hat=[]; 
%             end
%             
%             cdisp('Leave one out procedure started...');
%             blocks=[0 cumsum(obj.stem_data.dim)];
%             for j=1:obj.stem_data.nvar
%                 if sum(varset==j)>0
%                     residual{j}=nan(obj.stem_data.stem_varset.dim(j),obj.T);
%                     for i=1:obj.stem_data.stem_varset.dim(j)
%                         cdisp(['Site ',num2str(i),'/',num2str(obj.stem_data.stem_varset.dim(j)),' of ',obj.stem_data.stem_varset.name{j}]);
%                         obj.stem_par_initial=stem_par_alldata;
%                         chunk=obj.stem_data.Y(blocks(j)+i,:);
%                         obj.stem_data.set_Y_row(blocks(j)+i,chunk*NaN);
%                         obj.EM_estimate(EM_toll,EM_iter,'single');
%                         logL{j}(i)=obj.stem_EM_result.logL;
%                         chunk_full=obj.stem_EM_result.Y_hat(blocks(j)+i,:);
%                         residual{j}(i,:)=chunk-chunk_full;
%                         obj.stem_data.set_Y_row(blocks(j)+i,chunk);
%                     end
%                     rmse(j)=sum(sum(residual{j}(isnotnan(residual{j})).^2))/(sum(sum(isnotnan(residual{j}))));
%                     rmse(j)=sqrt(rmse(j));
%                 end
%             end
%             name=obj.stem_data.stem_varset.name;
%             for i=1:length(name)
%                 name{i}=['LOO Residuals - ',name{i}];
%             end
%             st_varsetstub=stem_varset(residual,name);
%             st_loo_residual=stem_residual(obj.stem_data,st_varsetstub);
%             st_loo_residual.kriging_Var_W_bar_hat=kriging_Var_W_bar_hat;
%             st_loo_residual.rmse=rmse;
%             st_loo_residual.logL=logL;
%             cdisp('Leave one out procedure ended');
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
        
        %set functions
        
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
                    disp('G and sigma_eta have been diagonalized');
                end
                obj.stem_par_initial=stem_par_initial;
            end
        end
        
        
       
    end
    
%     methods (Static)
%         
%         %far diventare il metodo non static?
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
    
end

