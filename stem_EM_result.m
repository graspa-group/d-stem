%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_EM_result < handle
    properties
        stem_par=[];                    %estimated stem_par object
        stem_par_all=[];                %list of stem_par objects of all EM iterations
        stem_kalmansmoother_result=[]   %stem_kalmansmoother_result object
        stem_residual=[];               %stem_residual object
        y_hat=[];                       %Estimated variable data with missing data filled
        y_hat_back=[];                  %Backtransformed estimated variable data with missing data filled
        res=[];                         %Model residuals
        res_back=[];                    %Backtransformed model residuals
        y_back=[];                      %Backtransformed original variable data
        E_wg_y1=[];
        E_wr_y1=[];
        Var_wg_y1=[];
        Var_wr_y1=[];
        computation_time=[];            %EM estimation time in seconds
        date_start=[];                  %estimation date
        machine=[];                     %estimation machine
        iterations=[];                  %number of EM iterations at convergence
        max_iterations=[];              %number of EM max iterations
        exit_toll=[];                   %EM exit toll
        logL=[];                        %log-likelihood value at convergence
        logL_all=[];                    %log-likelihood value at all iterations
        varcov=[];                      %parameter variance-covariance matrix at convergence
    end
    
    methods
        function obj = stem_EM_result()
            if nargin>0
                error('stem_EM_result constructor does not need input arguments');
            end
        end
    end
end