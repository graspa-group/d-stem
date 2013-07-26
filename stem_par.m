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

classdef stem_par
   
    properties
        %flags
        pixel_correlated=0;                 %[boolean]    (1x1) 1 if pixel sensing variables are correlated 0 otherwise        
        time_diagonal=0;                    %[boolean]    (1x1) 1: matrix G and sigma_eta are diagonal; 0: matrix G and sigma_eta are full
        clustering=0;                       %[boolean]    (1x1) 1: the property X_z of stem_data is estimated for clustering; 0: X_z is fixed and not estimated
        
        %fixed parameters
        q=[];                               %[integer]    (1x1) number of point level variables
        p=[];                               %[integer]    (1x1) dimension of the latent temporal state
        k=[];                               %[integer]    (1x1) number of coregionalization components for the point level variables
        n_beta=[];                          %[integer]    (1x1)
        correlation_type='exponential';     %[string]     (1x1) spatial correlation function type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2  
        theta_clustering=0;                 %[double >=0] (1x1) parameter of the spatial correlation fo clustering
        
        %estimated parameters
        beta=[];                            %[double]     (N_bx1) beta parameters
        alpha_bp=[];                        %[double]     (2qx1) alpha pixel/point parameters
        alpha_p=[];                         %[double]     (qx1) alpha point parameters;
        theta_b=[];                         %[double >0]  (q|1x1) coregionalization theta parameter for the pixel sensing variables
        theta_p=[];                         %[double >0]  (Kx1) coregionalization theta parameters for the point level variables
        v_b=[];                             %[double]     (qxq) correlation matrix for the pixel sensing variables
        v_p=[];                             %[double]     (qxqxK) correlation matrices for the point level variables
        G=[];                               %[double]     (pxp) transition matrix of the temporal component
        sigma_eta=[];                       %[double]     (pxp) error variance-covariance matrix of the temporal component
        sigma_eps=[];                       %[double]     (qxq|2qx2q) measurement-error variance-covariance matrix
    end
    
    methods
        
        function obj = stem_par(stem_data,correlation_type,pixel_correlated,time_diagonal,clustering,theta_clustering)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %[stem_data]                 - [stem_data object] (1x1)
            %<correlation_type>          - [string] (default:'exponential') (1x1) spatial correlation function type 
            %<pixel_correlated>          - [boolean] (default: 0) (1x1) 1: pixel variables are correlated; 0: otherwise 
            %<time_diagonal>             - [boolean] (default: 0) (1x1) 1: matrix G and sigma_eta are diagonal; 0: matrix G and sigma_eta are full
            %<clustering>                - [boolean] (dafault: 0) (1x1) 1: the property X_z of stem_data is estimated for clustering; 0: X_z is fixed and not estimated 
            %<theta_clustering>          - [double>=0] parameter of the spatial correlation for clustering
            %
            %OUTPUT
            %obj - [stem_par object] (1x1)
            
            if nargin<1
                error('The first input argument must be provided');
            end
            if not(isa(stem_data,'stem_data'))
                error('The first argument must be of class stem_data');
            end
            
            %q
            obj.q=length(stem_data.stem_varset_p.dim);
            
            %p
            tot=0;
            if not(isempty(stem_data.stem_varset_p.X_z))
                for i=1:length(stem_data.stem_varset_p.dim)
                    tot=tot+size(stem_data.stem_varset_p.X_z{i},2);
                end
            end
            if not(isempty(stem_data.stem_varset_b))
                if not(isempty(stem_data.stem_varset_b.X_z))
                    for i=1:length(stem_data.stem_varset_b.dim)
                        tot=tot+size(stem_data.stem_varset_b.X_z{i},2);
                    end
                end
            end
            obj.p=tot;
            if obj.p==0
                disp('No temporal component considered');
            end
            
            %n_beta
            tot=0;
            if not(isempty(stem_data.stem_varset_p.X_beta))
                for i=1:length(stem_data.stem_varset_p.dim)
                    tot=tot+size(stem_data.stem_varset_p.X_beta{i},2);
                end
            end
            if not(isempty(stem_data.stem_varset_b))
                if not(isempty(stem_data.stem_varset_b.X_beta))
                    for i=1:length(stem_data.stem_varset_b.dim)
                        tot=tot+size(stem_data.stem_varset_b.X_beta{i},2);
                    end
                end
            end
            obj.n_beta=tot;
            if obj.n_beta==0
                disp('No covariates considered');
            end            
            
            %k
            if isempty(stem_data.stem_varset_p.X_p)
                obj.k=0;
                disp('No geostatistical point level component considered');
            else
                obj.k=size(stem_data.stem_varset_p.X_p{1},4);
            end
            
            %correlation type
            if nargin<2
                disp('WARNING: Exponential correlation function is considered');
            else
                if not(isempty(correlation_type))
                    if sum(strcmp(correlation_type,{'exponential','matern32','matern52'}))==0
                        error('Only ''exponential'', ''matern32'' and ''matern52'' correlation functions are supported');
                    end
                    obj.correlation_type=correlation_type;
                else
                    disp('WARNING: the exponential correlation function is considered');
                end
            end
            
            if nargin>=3
                if not(isempty(pixel_correlated))
                    if isempty(stem_data.stem_varset_b)
                        disp('WARNING: Pixel data are not provided. The pixel_correlated input argument is ignored');
                    else
                        obj.pixel_correlated=pixel_correlated;
                    end
                end
            end
            
            if nargin>=4
                if not(isempty(time_diagonal))
                    if obj.p==0
                        disp('WARNING: p=0, the time_diagonal input argument is ignored');
                    else
                        obj.time_diagonal=time_diagonal;
                    end
                end
            end
            
            if nargin>=5
                if not(isempty(clustering))
                    if clustering==1
                        if obj.q>1
                            error('The clustering option is only available in the univariate case (q=1)');
                        end
                        if not(isempty(stem_data.stem_varset_b))
                            error('The clustering option is only available for point level data');
                        end
                        if isempty(stem_data.stem_varset_p.X_z)
                            error('X_z must be provided when the clustering option is enabled');
                        end
                    end
                    obj.clustering=clustering;
                end
            end
            
            if nargin>=6
                if not(isempty(theta_clustering))
                    if obj.clustering==0
                        disp('WARNING: the theta_clustering parameter is ignored');
                    else
                        obj.theta_clustering=theta_clustering;
                    end
                end
            end
            
            %matrix building
            
            if obj.n_beta>0
                obj.beta=zeros(obj.n_beta,1);
            end
            if not(isempty(stem_data.stem_varset_b))
                obj.alpha_bp=zeros(obj.q*2,1);
                obj.sigma_eps=zeros(obj.q*2);
                if obj.pixel_correlated
                    obj.theta_b=0;
                else
                    obj.theta_b=zeros(obj.q,1);
                end                
                obj.v_b=eye(obj.q);
            else
                obj.sigma_eps=zeros(obj.q);
            end
            
            if obj.k>0
                obj.alpha_p=zeros(obj.q,obj.k);
                obj.theta_p=zeros(obj.k,1);
                for i=1:obj.k
                    v_p(:,:,i)=eye(obj.q);
                end
                obj.v_p=v_p;
            end
            if obj.p>0
                obj.G=zeros(obj.p);
                obj.sigma_eta=zeros(obj.p);
            end
        end
        
        function all_par = vec(obj)
            %DESCRIPTION: vectorize the model parameters (only the estimated parameters with the exclusion of structural zeroes and repeated elements of symmetric matrices)
            %
            %INPUT
            %obj     - [stem_par object] (1x1) 
            %
            %OUTPUT
            %all_par - [double]          (Hx1)
                  
            all_par=[];
            all_par=[all_par; obj.beta];
            all_par=[all_par; diag(obj.sigma_eps)];

            all_par=[all_par; obj.alpha_bp];
            all_par=[all_par; obj.theta_b(:)];

            if obj.pixel_correlated
                all_par=[all_par; stem_misc.from_upper_triangular_to_vector(obj.v_b)];
            end

            all_par=[all_par; obj.alpha_p(:)];
            all_par=[all_par;obj.theta_p(:)];
            for i=1:obj.k
                all_par=[all_par; stem_misc.from_upper_triangular_to_vector(obj.v_p(:,:,i))];
            end
            
            if obj.p>0
                if not(obj.time_diagonal)
                    all_par=[all_par; obj.G(:)];
                    all_par=[all_par; stem_misc.triuv(obj.sigma_eta)];
                else
                    all_par=[all_par; diag(obj.G)];
                    all_par=[all_par; diag(obj.sigma_eta)];
                end
            end
        end
        
        function print(obj)
            %DESCRIPTION: print the stem_par object
            %
            %INPUT
            %obj     - [stem_par object] (1x1) 
            %
            %OUTPUT
            %none: the stem_par object is printed in the command window        
            disp('****************');
            disp('PARAMETER VALUES');
            disp('****************');
            if not(isempty(obj.alpha_bp))
                disp(['alpha_bp: ',num2str(obj.alpha_bp')]);
            end
            if not(isempty(obj.alpha_p))
            disp(['alpha_p: ',num2str(obj.alpha_p(:)')]);
            end
            if not(isempty(obj.theta_b))
            disp(['theta_b: ',num2str(obj.theta_b(:)')]);
            end
            if not(isempty(obj.theta_p))
                disp(['theta_p: ',num2str(obj.theta_p(:)')]);
            end
            if not(isempty(obj.v_b))
                if (obj.pixel_correlated)&&(obj.q>1)
                    disp(['v_b: ']);
                    disp(obj.v_b);
                end
            end
            if obj.q>1
                for i=1:obj.k
                    disp(['v_p',num2str(i),': ']);
                    disp(obj.v_p(:,:,i));
                end
            end
            if not(isempty(obj.beta))
            disp(['beta: ',num2str(obj.beta')]);
            end
            disp(['sigma eps: ',num2str(diag(obj.sigma_eps)')]);
            if obj.p>0
                if obj.time_diagonal
                    disp(['diag_G: ',num2str(diag(obj.G)')]);
                    disp(['sigma_eta: ',num2str(diag(obj.sigma_eta)')]);
                else
                    disp(['G: ']);
                    disp(obj.G);
                    disp(['sigma_eta: ']);
                    disp(obj.sigma_eta);
                end
            end
        end

        %Class set methods
        function obj = set.q(obj,q)
            obj.q=q;
        end
        
        function obj = set.k(obj,k)
            if k<0
                error('k must be >=0');
            end
            obj.k=k;
        end
        
        function obj = set.beta(obj,beta)
            obj.beta=beta;
        end
        
        function obj = set.alpha_bp(obj,alpha_bp)
            if length(alpha_bp)~=obj.q*2
                error(['The length of alpha_bp must be equal to ',num2str(2*obj.q)]);
            end
            obj.alpha_bp=alpha_bp;
        end
        
        function obj = set.alpha_p(obj,alpha_p)
            if not(size(alpha_p,1)==obj.q && size(alpha_p,2)==obj.k)
                error(['alpha_p must be ',num2str(obj.q),'x',num2str(obj.k)]);
            end
            obj.alpha_p=alpha_p;
        end        
        
        function obj = set.theta_b(obj,theta_b)
            if not(obj.pixel_correlated) && not(length(theta_b)==obj.q)
                error(['The length of theta_b must be equal to ',num2str(obj.q)]);
            end
            if obj.pixel_correlated && not(length(theta_b)==1)
                error('theta_b must be a scalar');
            end
            if sum(theta_b<0)>0
                error('The element of theta_b cannot be negative');
            end
            obj.theta_b=theta_b;
        end
        
        function obj = set.theta_p(obj,theta_p)
            if not(length(theta_p)==obj.k)
                error(['theta_p must be 1x',num2str(obj.k)]);
            end
            if sum(theta_p<0)>0
                error('The element of theta_p cannot be negative');
            end
            obj.theta_p=theta_p;
        end

        function obj = set.v_b(obj,v_b)
            if not(size(v_b,1)==obj.q)||(size(v_b,2)~=obj.q)
                error('v_b must be qxq');
            end
            if not(obj.pixel_correlated)
                temp=v_b-eye(size(v_b,1));
                if sum(temp(:))>0
                    error('v_b must be the identity matrix since the pixel variables are uncorrelated');
                end
            else
                if not(sum(diag(v_b-eye(size(v_b,1))))==0)
                    error('The diagonal elements of v_b must be 1');
                end                
            end
            if min(eig(v_b))<0
                error('v_b must be positive definited');
            end
            obj.v_b=v_b;
        end
        
        function obj = set.v_p(obj,v_p)
            if not(size(v_p,1)==obj.q)||not(size(v_p,2)==obj.q)||not(size(v_p,3)==obj.k)
                error(['v_p must be ',num2str(obj.q),'x',num2str(obj.q),'x',num2str(obj.k)]);
            end
            for i=1:size(v_p,3)
                if not(sum(diag(v_p(:,:,i)-eye(size(v_p,1))))==0)
                    error('The diagonal elements of each v_p(:,:,i) matrix must be 1');
                end
                if min(eig(v_p(:,:,i)))<0
                    error('Each v_p(:,:,i) matrix must be positive definited');
                end
            end
            obj.v_p=v_p;
        end        
        
        function obj = set.G(obj,G)
            if not(size(G,1)==obj.p) || not(size(G,2)==obj.p)
                error(['G must be a ',num2str(obj.p),'x',num2str(obj.p),' matrix']);
            end
            obj.G=G;
        end
        
        function obj = set.sigma_eta(obj,sigma_eta)
            if not(size(sigma_eta,1)==obj.p) || not(size(sigma_eta,2)==obj.p)
                error(['sigma_eta must be a ',num2str(obj.p),'x',num2str(obj.p),' matrix']);
            end
            obj.sigma_eta=sigma_eta;            
        end
        
        function obj = set.sigma_eps(obj,sigma_eps)
            if not(size(sigma_eps,1)==size(sigma_eps,2))
                error('sigma_eps must be square');
            end        
            obj.sigma_eps=sigma_eps;
        end
        
        function obj = set.theta_clustering(obj,theta_clustering)
            if theta_clustering<0
                error('theta_clustering must be >= 0')
            end
            obj.theta_clustering=theta_clustering;
        end        
        
    end
end