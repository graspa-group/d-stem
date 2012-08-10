classdef stem_par_constraint
    properties
        sigma_eps_diag=1;
        sigma_eta_diag=0;
        no_direct_components=0;
        K_from_data=0;
    end
    
    methods
        function obj = stem_par_constraint()
            if nargin>0
                error('The constructor does not need input arguments');
            end
        end
        
        function obj = set.sigma_eps_diag(obj,sigma_eps_diag)
            if sum(sigma_eps_diag==[0 1])==0
                error('sigma_eps_diag must be either 0 or 1');
            end
            obj.sigma_eps_diag=sigma_eps_diag;
        end
        
        function obj = set.sigma_eta_diag(obj,sigma_eta_diag)
            if sum(sigma_eta_diag==[0 1])==0
                error('sigma_eta_diag must be either 0 or 1');
            end
            obj.sigma_eta_diag=sigma_eta_diag;
        end
        
        function obj = set.no_direct_components(obj,no_direct_components)
            if sum(no_direct_components==[0 1])==0
                error('no_direct_components must be either 0 or 1');
            end
            obj.no_direct_components=no_direct_components;
        end
        
        function obj = set.K_from_data(obj,K_from_data)
            if sum(K_from_data==[0 1])==0
                error('K_from_data must be either 0 or 1');
            end
            obj.K_from_data=K_from_data;
        end
        
    end
end