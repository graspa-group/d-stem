%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_kalmanfilter_result < handle
    properties
        zk_f = [];
        zk_u = [];
        Pk_f = [];
        Pk_u = [];
        J    = [];
        logL = [];
    end
    
    methods
        function obj = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J,logL)
            obj.zk_f=zk_f;
            obj.zk_u=zk_u;
            obj.Pk_f=Pk_f;
            obj.Pk_u=Pk_u;
            obj.J=J;
            obj.logL=logL;
        end
    end
end