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
        %estimated parameters
        beta=[];            %[double]       (N_bx1) beta parameters
        alpha_bp=[];        %[double]       (2qx1) alpha pixel/point parameters
        alpha_p=[];         %[double]       (qxk) alpha point parameters;
        theta_b=[];         %[double >0]    (qx1|1x1) coregionalization theta parameter for the pixel variables
        theta_p=[];         %[double >0]    (kx1|kxd) coregionalization theta parameters for the point level variables
        v_b=[];             %[double]       (qxq) correlation matrix for the pixel sensing variables
        v_p=[];             %[double]       (qxqxk) correlation matrices for the point level variables
        G=[];               %[double]       (pxp) transition matrix of the temporal component
        lambda=[];          %[double]       (1x1) lambda parameter for the case of irregular time steps
        sigma_eta=[];       %[double]       (pxp) error variance-covariance matrix of the temporal component
        sigma_eps=[];       %[double]       (qxq|2qx2q|qxqxT) measurement-error variance-covariance matrix
        alpha_z=[];         %[double]       (px1) alpha parameters related to the z latent variable when model_name is 'HDGM' or 'f-HDGM'
        theta_z=[];         %[double]       (1x1|1x2|1xp|2xp) coregionalization theta parameter related to the z latent variable when model_name is 'HDGM' or 'f-HDGM'
        v_z=[];             %[double]       (pxp) coregionalization matrix related to the z latent variable when model_name is 'HDGM' or 'f-HDGM'
        
        stem_fda=[];        %[stem_fda obj] (1x1) object of class stem_fda
    end
    
    properties (SetAccess=private)
        %model type
        stem_modeltype=[];              %[stem_modeltype obj] (1x1) object of class stem_modeltype
        %par constraints
        stem_par_constraints=[];        %[stem_par_contraints obj] (1x1) object of class stem_par_contraints 
        
        %fixed parameters
        p=[];                           %[integer>0]    (1x1) dimension of the latent temporal variable z
        q=[];                           %[integer>0]    (1x1) number of point level variables
        k=[];                           %[integer>0]    (1x1) number of coregionalization components for the point level variables
        n_beta=[];                      %[integer>0]    (1x1) total number of loading coefficients in X_beta
        d=[];                           %[integer>0]    (1x1)
        correlation_type='exponential'; %[string]       (1x1) spatial correlation function type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2; 'expsphere': anisotropic correlation function on the sphere  
    end
    
    methods
        
        function obj = stem_par(obj_stem_data,correlation_type,obj_stem_fda,obj_stem_par_constraints)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %obj_stem_data                  - [stem_data object]            (1x1) an object of class stem_data
            %<correlation_type>             - [string]                      (1x1) spatial correlation function type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2; 'expsphere': anisotropic correlation function on the sphere
            %<obj_stem_fda>                 - [stem_fda object]             (1x1) an object of class stem_fda
            %<obj_stem_par_contraints>      - [stem_par_contraints object]	(1x1) an object of class stem_par_constraints
            %
            %OUTPUT
            %obj                        - [stem_par object] (1x1)
            
            if nargin<1
                error('stem_data must be provided');
            end
            if not(isa(obj_stem_data,'stem_data'))
                error('obj_stem_data must be of class stem_data');
            end
            if nargin>1
                obj.correlation_type=correlation_type;
            end
            if nargin>2
                if not(isa(obj_stem_fda,'stem_fda'))
                    error('obj_stem_fda must be of class stem_fda');
                end
                obj.stem_fda=obj_stem_fda;
            end
            if nargin>3
                if not(isa(obj_stem_par_constraints,'stem_par_constraints'))
                    error('obj_stem_par_constraints must be of class stem_par_constraints');
                end
                obj.stem_par_constraints=obj_stem_par_constraints;
            else
                obj.stem_par_constraints=stem_par_constraints();
            end
            
            %stem_modeltype
            obj.stem_modeltype=obj_stem_data.stem_modeltype;

            %d
            obj.d=size(obj_stem_data.stem_gridlist_p.grid{1}.coordinate,2);
            
            %q
            obj.q=length(obj_stem_data.stem_varset_p.dim);
            
            %p
            if not(obj.stem_modeltype.is({'HDGM','f-HDGM'}))
                tot=0;
                if not(isempty(obj_stem_data.stem_varset_p.X_z))
                    for i=1:length(obj_stem_data.stem_varset_p.dim)
                        if sum(abs(obj_stem_data.stem_varset_p.X_z{i}(:)))>0
                            tot=tot+size(obj_stem_data.stem_varset_p.X_z{i},2);
                        end
                    end
                end
                if not(isempty(obj_stem_data.stem_varset_b))
                    if not(isempty(obj_stem_data.stem_varset_b.X_z))
                        for i=1:length(obj_stem_data.stem_varset_b.dim)
                            if sum(abs(obj_stem_data.stem_varset_b.X_z{i}(:)))>0
                                tot=tot+size(obj_stem_data.stem_varset_b.X_z{i},2);
                            end
                        end
                    end
                end
                obj.p=tot;
            else
                if obj.stem_modeltype.is('HDGM')
                    if not(isempty(obj_stem_data.stem_varset_p.X_z))
                        if size(obj_stem_data.stem_varset_p.X_z{1},2)>1
                            obj.p=size(obj_stem_data.stem_varset_p.X_z{1},2);
                        else
                            obj.p=obj.q;
                        end
                    end
                else
                    if nargin<3
                        error('obj_stem_fda must be provided since model_type is f-HDGM');
                    end
                    obj.p=getnbasis(obj.stem_fda.spline_basis);
                end
            end
            
            %n_beta
            if not(obj.stem_modeltype.is('f-HDGM'))
                tot=0;
                if not(isempty(obj_stem_data.stem_varset_p.X_beta))
                    for i=1:length(obj_stem_data.stem_varset_p.dim)
                        if sum(abs(obj_stem_data.stem_varset_p.X_beta{i}(:)))>0
                            tot=tot+size(obj_stem_data.stem_varset_p.X_beta{i},2);
                        end
                    end
                end
                if not(isempty(obj_stem_data.stem_varset_b))
                    if not(isempty(obj_stem_data.stem_varset_b.X_beta))
                        for i=1:length(obj_stem_data.stem_varset_b.dim)
                            if sum(abs(obj_stem_data.stem_varset_b.X_beta{i}(:)))>0
                                tot=tot+size(obj_stem_data.stem_varset_b.X_beta{i},2);
                            end
                        end
                    end
                end
                obj.n_beta=tot;
            else
                size_ref=0;
                if not(isempty(obj_stem_data.stem_varset_p.X_beta))
                    size_ref=size(obj_stem_data.stem_varset_p.X_beta{1},2);
                    for i=1:length(obj_stem_data.stem_varset_p.dim)
                        if sum(abs(obj_stem_data.stem_varset_p.X_beta{i}(:)))>0
                            if not(size(obj_stem_data.stem_varset_p.X_beta{i},2)==size_ref)
                                error('When shared_beta=1, all the X_beta cells must have the same covariates');
                            end
                        end
                    end
                end
                if not(isempty(obj_stem_data.stem_varset_b))
                    if not(isempty(obj_stem_data.stem_varset_b.X_beta))
                        for i=1:length(obj_stem_data.stem_varset_b.dim)
                            if sum(abs(obj_stem_data.stem_varset_b.X_beta{i}(:)))>0
                                if not(size(obj_stem_data.stem_varset_b.X_beta{i},2)==size_ref)
                                    error('When model_name is ''f-HDGM'', all the X_beta cells must have the same covariates');
                                end
                            end
                        end
                    end
                end
                obj.n_beta=size_ref;
            end

            %k
            if isempty(obj_stem_data.stem_varset_p.X_p)
                obj.k=0;
            else
                obj.k=size(obj_stem_data.stem_varset_p.X_p{1},4);
            end
            
            %parameter vector and matrix building
            
            if obj.n_beta>0
                obj.beta=zeros(obj.n_beta,1);
            end
            
            
            if not(obj.stem_modeltype.is('MBC')&&strcmpi(obj.stem_modeltype.clustering_type,'Dynamic'))
                if not(obj.stem_modeltype.is('f-HDGM'))
                    if not(isempty(obj_stem_data.stem_varset_b))
                        obj.alpha_bp=zeros(obj.q*2,1);
                        obj.sigma_eps=zeros(obj.q*2);
                        if obj.stem_par_constraints.pixel_correlated
                            if not(strcmp(obj.correlation_type,'expsphere'))
                                obj.theta_b=0;
                            else
                                obj.theta_b=zeros(2,1);
                            end
                        else
                            if not(strcmp(obj.correlation_type,'expsphere'))
                                obj.theta_b=zeros(1,obj.q);
                            else
                                obj.theta_b=zeros(2,obj.q);
                            end
                        end
                        obj.v_b=eye(obj.q);
                    else
                        obj.sigma_eps=zeros(obj.q);
                    end
                else
                    obj.sigma_eps=0;
                end
            else
                obj.sigma_eps=zeros(obj.p);
            end
            
            if obj.k>0
                obj.alpha_p=zeros(obj.q,obj.k);
                if not(obj.stem_modeltype.is('Emulator'))
                    if not(strcmp(obj.correlation_type,'expsphere'))
                        obj.theta_p=zeros(1,obj.k);
                    else
                        obj.theta_p=zeros(2,obj.k);
                    end
                else
                    obj.theta_p=zeros(obj.k,obj.d);
                end
                vp=zeros(obj.q,obj.q,obj.k);
                for i=1:obj.k
                    vp(:,:,i)=eye(obj.q);
                end
                obj.v_p=vp;
            end
            if obj.p>0
                obj.G=zeros(obj.p);
                if not(obj.stem_modeltype.is({'HDGM','f-HDGM'}))
                    obj.sigma_eta=zeros(obj.p);
                else
                    obj.alpha_z=zeros(obj.p,1);
                    if not(strcmp(correlation_type,'expsphere'))
                        if obj.stem_modeltype.is('f-HDGM')
                            obj.theta_z=zeros(1,obj.p);
                        else
                            obj.theta_z=0;
                        end
                    else
                        if obj.stem_modeltype.is('f-HDGM')
                            obj.theta_z=zeros(2,obj.p);
                        else
                            obj.theta_z=[0 0]';
                        end
                    end
                    obj.v_z=eye(obj.p);
                end
                if obj_stem_data.stem_datestamp.irregular==1
                    obj.lambda=1;
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
            stem_misc.compare(obj.beta,beta,stem_misc.varname(beta));
            obj.beta=beta;
        end
        
        function obj = set.alpha_bp(obj,alpha_bp)
            stem_misc.compare(obj.alpha_bp,alpha_bp,stem_misc.varname(alpha_bp));
            obj.alpha_bp=alpha_bp;
        end
        
        function obj = set.alpha_p(obj,alpha_p)
            stem_misc.compare(obj.alpha_p,alpha_p,stem_misc.varname(alpha_p));
            obj.alpha_p=alpha_p;
        end      
        
        function obj = set.alpha_z(obj,alpha_z)
            stem_misc.compare(obj.alpha_z,alpha_z,stem_misc.varname(alpha_z));
            obj.alpha_z=alpha_z;
        end 
        
        function obj = set.theta_b(obj,theta_b)
            if sum(theta_b<0)>0
                error('The element of theta_b cannot be negative');
            end
            stem_misc.compare(obj.theta_b,theta_b,stem_misc.varname(theta_b));
            obj.theta_b=theta_b;
        end
        
        function obj = set.theta_p(obj,theta_p)
            if sum(theta_p(:)<0)>0
                error('The element of theta_p cannot be negative');
            end
            stem_misc.compare(obj.theta_p,theta_p,stem_misc.varname(theta_p));
            obj.theta_p=theta_p;
        end
        
        function obj = set.theta_z(obj,theta_z)
            stem_misc.compare(obj.theta_z,theta_z,stem_misc.varname(theta_z));
            obj.theta_z=theta_z;
        end        
        
        function obj = set.v_b(obj,v_b)
            if not(obj.stem_par_constraints.pixel_correlated)
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
            stem_misc.compare(obj.v_b,v_b,stem_misc.varname(v_b));
            obj.v_b=v_b;
        end
        
        function obj = set.v_p(obj,v_p)
            for i=1:size(v_p,3)
                if not(sum(diag(v_p(:,:,i)-eye(size(v_p,1))))==0)
                    error('The diagonal elements of each v_p(:,:,i) matrix must be 1');
                end
                if min(eig(v_p(:,:,i)))<0
                    error('Each v_p(:,:,i) matrix must be positive definited');
                end
            end
            stem_misc.compare(obj.v_p,v_p,stem_misc.varname(v_p));
            obj.v_p=v_p;
        end   
        
        function obj = set.v_z(obj,v_z)
            if not(sum(diag(v_z-eye(size(v_z,1))))==0)
                error('The diagonal elements of v_z must be 1');
            end
            if min(eig(v_z))<0
                error('v_z must be positive definited');
            end
            stem_misc.compare(obj.v_z,v_z,stem_misc.varname(v_z));
            obj.v_z=v_z;
        end        
        
        function obj = set.G(obj,G)
            stem_misc.compare(obj.G,G,stem_misc.varname(G));
            obj.G=G;
        end
        
        function obj = set.sigma_eta(obj,sigma_eta)
            stem_misc.compare(obj.sigma_eta,sigma_eta,stem_misc.varname(sigma_eta));
            obj.sigma_eta=sigma_eta;            
        end
        
        function obj = set.sigma_eps(obj,sigma_eps)
            stem_misc.compare(obj.sigma_eps,sigma_eps,stem_misc.varname(sigma_eps));
            obj.sigma_eps=sigma_eps;
        end
        
        function obj = set.stem_fda(obj,stem_fda)
            if not(isa(stem_fda,'stem_fda'))
                error('stem_fda must be of class stem_fda');
            end
            obj.stem_fda=stem_fda;
        end
    end
end