%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-STEM - Distributed Space Time Expecation Maximization      %
%                                                              %
% Author: Francesco Finazzi                                    %
% E-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo - Dept. of Engineering    %
% Author website: http://www.unibg.it/pers/?francesco.finazzi  %
% Code website: https://code.google.com/p/d-stem/              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



classdef stem_EM_result < handle
    
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %T   - number of temporal steps
    %H   - total number of scalar model parameters
    %I   - number of EM iterations at convergence 
    
    properties
        stem_par=[];                    %[stem_par object]  (1x1) estimated model parameters after the last EM itaration
        stem_par_all=[];                %[stem_par object]  (Ix1) estimated model parameters at all the EM iterations
        stem_kalmansmoother_result=[]   %[stem_kalmansmoother_result object]    (1x1) stem_kalmansmoother_result object
   
        y_hat=[];                       %[double]           (NxT) Estimated variable data with missing data filled
        y_hat_back=[];                  %[double]           (NxT) Backtransformed estimated variable data with missing data filled
        y_back=[];                      %[double]           (NxT) Backtransformed original variable data
        res=[];                         %[double]           (NxT) Model residuals
        res_back=[];                    %[double]           (NxT) Backtransformed model residuals
        E_wp_y1=[];                     %[double]           (N_pxTxK) E[wp|Y(1)]
        E_wb_y1=[];                     %[double]           (N_pxTxK) E[wp_k|Y(1)]
        diag_Var_wp_y1=[];              %[double]           (N_pxTxK) diagonals of Var[wp|Y(1)]
        diag_Var_wb_y1=[];              %[double]           (N_bxT)   diagonals of Var[wb|Y(1)]
        
        logL=[];                        %[double]           (1x1) observed-data log-likelihood at convergence
        logL_all=[];                    %[double]           (Ix1) observed-data log-likelihood at all EM iterations
        varcov=[];                      %[double]           (HxH) parameter variance-covariance matrix at convergence
        
        iterations=[];                  %[integer >0]       (1x1) number of EM iterations at convergence
        max_iterations=[];              %[integer >0]       (1x1) maximum number of EM iterations
        exit_toll=[];                   %[double] >0]       (1x1) exit tolerance
        
        computation_time=[];            %[double]           (1x1) EM estimation time in seconds
        date_start=[];                  %[double]           (1x1) date of the start of the EM estimation
        machine=[];                     %[string]           (1x1) operative system of the machine
    end
    
    methods
        function obj = stem_EM_result()
            if nargin>0
                error('stem_EM_result constructor does not need input arguments');
            end
        end
    end
end