%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_par
    %class for stem parameters
   
    properties
        %flag
        remote_correlated=[];               %[1 x 1 boolean] 1 if remote sensing variables are correlated 0 otherwise        
        
        %fixed parameters
        q=[];                               %[1 x 1 integer] number of ground level variables
        p=[];                               %[1 x 1 integer] dimension of the latent temporal state
        k=[];                               %[1 x 1 integer] number of coregionalization components for the ground level variables
        n_beta=[];                          %[1 x 1 integer]
        correlation_type='';                %[string] correlation function type  
        %constraints
        time_diagonal=0;                    %[1 x 1 boolean]
        
        %estimated parameters
        beta=[];                            %[dim(n_beta) x 1 double] beta parameters
        alpha_rg=[];                        %[2q x 1 double] alpha remote/ground parameters
        alpha_g=[];                         %[q x 1 double] alpha ground parameters;
        theta_r=[];                         %[1 x 1 or q x 1 double] coregionalization theta parameter for the remote sensing variables
        theta_g=[];                         %[k x 1 double] coregionalization theta parameters for the ground level variables
        v_r=[];                             %[q x q double] correlation matrix for the remote sensing variables
        v_g=[];                             %[q x q x k double] correlation matrices for the ground level variables
        G=[];                               %[p x p double] transition matrix of the temporal component
        sigma_eta=[];                       %[p x p double] error variance-covariance matrix of the temporal component
        sigma_eps=[];                       %[q x q | 2q x 2q] measurement-error variance-covariance matrix
    end
    
    methods
        
        function obj = stem_par(stem_data,correlation_type,remote_correlated,time_diagonal)
            %costructor
            if nargin<1
                error('The first input argument must be provided');
            end
            if not(isa(stem_data,'stem_data'))
                error('The first argument must be of class stem_data');
            end
            
            %q
            obj.q=length(stem_data.stem_varset_g.dim);
            
            %p
            tot=0;
            if not(isempty(stem_data.stem_varset_g.X_time))
                for i=1:length(stem_data.stem_varset_g.dim)
                    tot=tot+size(stem_data.stem_varset_g.X_time{i},2);
                end
            end
            if not(isempty(stem_data.stem_varset_r))
                if not(isempty(stem_data.stem_varset_r.X_time))
                    for i=1:length(stem_data.stem_varset_r.dim)
                        tot=tot+size(stem_data.stem_varset_r.X_time{i},2);
                    end
                end
            end
            obj.p=tot;
            if obj.p==0
                disp('No temporal component considered');
            end
            
            %n_beta
            tot=0;
            if not(isempty(stem_data.stem_varset_g.X_beta))
                for i=1:length(stem_data.stem_varset_g.dim)
                    tot=tot+size(stem_data.stem_varset_g.X_beta{i},2);
                end
            end
            
            if not(isempty(stem_data.stem_varset_r))
                if not(isempty(stem_data.stem_varset_r.X_beta))
                    for i=1:length(stem_data.stem_varset_r.dim)
                        tot=tot+size(stem_data.stem_varset_r.X_beta{i},2);
                    end
                end
            end
            obj.n_beta=tot;
            if obj.n_beta==0
                disp('No covariates considered');
            end            
            
            %k
            if isempty(stem_data.stem_varset_g.X_g)
                obj.k=0;
                disp('No geostatistical ground level component considered');
            else
                obj.k=size(stem_data.stem_varset_g.X_g{1},4);
            end
            
            %correlation type
            if nargin<2
                obj.correlation_type='exponential';
                disp('WARNING: Exponential correlation function is considered');
            else
                if not(isempty(correlation_type))
                    if sum(strcmp(correlation_type,{'exponential','matern'}))==0
                        error('Only ''exponential'' and ''matern'' correlation functions are supported');
                    end
                    obj.correlation_type=correlation_type;
                else
                    obj.correlation_type='exponential';
                    disp('WARNING: Exponential correlation function is considered');
                end
            end
            
            if nargin>=3
                if isempty(stem_data.stem_varset_r)
                    disp('WARNING: Remote sensing data are not provided. The remote_correlated input argument is ignored');
                else
                    obj.remote_correlated=remote_correlated;
                end
            end
            
            if nargin>=4
                if obj.p==0
                    disp('WARNING: p=0, the time_diagonal input argument is ignored');
                else
                    obj.time_diagonal=time_diagonal;
                end
            end
            
            
            %matrix building
            
            if obj.n_beta>0
                obj.beta=zeros(obj.n_beta,1);
            end
            if not(isempty(stem_data.stem_varset_r))
                obj.alpha_rg=zeros(obj.q*2,1);
                obj.sigma_eps=zeros(obj.q*2);
                if strcmp(obj.correlation_type,'exponential') && obj.remote_correlated
                    obj.theta_r=0;
                end
                if strcmp(obj.correlation_type,'matern') && obj.remote_correlated
                    obj.theta_r=zeros(1,2);
                end
                if strcmp(obj.correlation_type,'exponential') && not(obj.remote_correlated)
                    obj.theta_r=zeros(obj.q,1);
                end                
                if strcmp(obj.correlation_type,'matern') && not(obj.remote_correlated)
                    obj.theta_r=zeros(obj.q,2);
                end 
                obj.v_r=eye(obj.q);
            else
                obj.sigma_eps=zeros(obj.q);
            end
            
            if obj.k>0
                obj.alpha_g=zeros(obj.q,obj.k);
                if strcmp(obj.correlation_type,'exponential')
                    obj.theta_g=zeros(obj.k,1);
                else
                    obj.theta_g=zeros(obj.k,2);
                end
                
                for i=1:obj.k
                    v_g(:,:,i)=eye(obj.q);
                end
                obj.v_g=v_g;
            end
            if obj.p>0
                obj.G=zeros(obj.p);
                obj.sigma_eta=zeros(obj.p);
            end
        end
        
        function all_par = vec(obj)
            all_par=[];
            all_par=[all_par; obj.beta];
            all_par=[all_par; diag(obj.sigma_eps)];

            all_par=[all_par; obj.alpha_rg];
            if strcmp(obj.correlation_type,'exponential')
                all_par=[all_par; obj.theta_r];
            else
                all_par=[all_par; obj.theta_r(:)];
            end
            if obj.remote_correlated
                all_par=[all_par; stem_par.from_upper_triangular_to_vector(obj.v_r)];
            end

            all_par=[all_par; obj.alpha_g(:)];
            if strcmp(obj.correlation_type,'exponential')
                all_par=[all_par; obj.theta_g];
            else
                all_par=[all_par;obj.theta_g(:)];
            end
            for i=1:obj.k
                all_par=[all_par; stem_par.from_upper_triangular_to_vector(obj.v_g(:,:,i))];
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
            disp('****************');
            disp('PARAMETER VALUES');
            disp('****************');
            disp(['alpha_rg: ',num2str(obj.alpha_rg')]);
            disp(['alpha_g: ',num2str(obj.alpha_g)]);
            disp(['theta_r: ',num2str(obj.theta_r')]);
            disp(['theta_g: ',num2str(obj.theta_g')]);
            if not(isempty(obj.v_r))
                if (obj.remote_correlated)&&(obj.q>1)
                    disp(['v_r: ']);
                    disp(obj.v_r);
                end
            end
            if obj.q>1
                for i=1:obj.k
                    disp(['v_g',num2str(i),': ']);
                    disp(obj.v_g(:,:,i));
                end
            end
            disp(['beta: ',num2str(obj.beta')]);
            disp(['sigma eps: ',num2str(diag(obj.sigma_eps)')]);
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
        
        function obj = set.alpha_rg(obj,alpha_rg)
            if length(alpha_rg)~=obj.q*2
                error(['The length of alpha_rg must be equal to ',num2str(2*obj.q)]);
            end
            obj.alpha_rg=alpha_rg;
        end
        
        function obj = set.alpha_g(obj,alpha_g)
            if not(size(alpha_g,1)==obj.q && size(alpha_g,2)==obj.k)
                error(['alpha_g must be ',num2str(obj.q),'x',num2str(obj.k)]);
            end
            obj.alpha_g=alpha_g;
        end        
        
        function obj = set.theta_r(obj,theta_r)
            if strcmp(obj.correlation_type,'exponential')
                if not(obj.remote_correlated) && not(length(theta_r)==obj.q)
                    error(['The length of theta_r must be equal to ',num2str(obj.q)]);
                end
                if obj.remote_correlated && not(length(theta_r)==1)
                    error('theta_r must be a scalar');
                end
                if sum(theta_r<0)>0
                    error('The element of theta_r cannot be negative');
                end
            end
            if strcmp(obj.correlation_type,'matern')
                if not(obj.remote_correlated) && not(size(theta_r,2)==2) && not(size(theta_r,1)==obj.q*2)
                    error(['The size of theta_r must be ',num2str(2*obj.q),' x 2']);
                end
                if obj.remote_correlated && not(size(theta_r,2)==2) && not(size(theta_r,1)==1)
                    error('The size of theta_r must be 1 x 2');
                end                
            end
            obj.theta_r=theta_r; 
        end
        
        function obj = set.theta_g(obj,theta_g)
            if strcmp(obj.correlation_type,'exponential')
                if not(length(theta_g)==obj.k)
                    error(['The length of theta_g must be equal to ',num2str(obj.k)]);
                end
                if sum(theta_g<0)>0
                    error('The element of theta_g cannot be negative');
                end
            end
            if strcmp(obj.correlation_type,'matern')
                if not(size(theta_g,2)==2) && not(size(theta_g,1)==obj.k)
                    error(['The size of theta_g must be ',num2str(obj.k),' x 2']);
                end
            end
            obj.theta_g=theta_g; 
        end        

        function obj = set.v_r(obj,v_r)
            if not(size(v_r,1)==obj.q)||(size(v_r,2)~=obj.q)
                error('v_r must be qxq');
            end
            if not(obj.remote_correlated)
                temp=v_r-eye(size(v_r,1));
                if sum(temp(:))>0
                    error('v_r must be the identity matrix since the remote variables are uncorrelated');
                end
            else
                if not(sum(diag(v_r-eye(size(v_r,1))))==0)
                    error('The diagonal elements of v_r must be 1');
                end                
            end
            if min(eig(v_r))<0
                error('v_r must be positive definited');
            end
            obj.v_r=v_r;
        end
        
        function obj = set.v_g(obj,v_g)
            if not(size(v_g,1)==obj.q)||not(size(v_g,2)==obj.q)||not(size(v_g,3)==obj.k)
                error(['v_g must be ',num2str(obj.q),'x',num2str(obj.q),'x',num2str(obj.k)]);
            end
            for i=1:size(v_g,3)
                if not(sum(diag(v_g(:,:,i)-eye(size(v_g,1))))==0)
                    error('The diagonal elements of each v_g(:,:,i) matrix must be 1');
                end
                if min(eig(v_g(:,:,i)))<0
                    error('Each v_g(:,:,i) matrix must be positive definited');
                end
            end
            obj.v_g=v_g;
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
    end
    
    methods (Static)
        function vec = from_upper_triangular_to_vector(mat)
            d=size(mat,1);
            vec=zeros(d*(d-1)/2,1);
            
            counter=1;
            for i=1:d-1
                for j=i+1:d
                    vec(counter)=mat(i,j);
                    counter=counter+1;
                end
            end
            
        end
    end
end