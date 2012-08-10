%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef EM < handle
    %abstract EM class
    
    properties (Abstract)
        stem_EM_options;
    end
    
    methods (Abstract)
        estimate(varargin)
        varargout = E_step(varargin)
        varargout = M_step(varargin)
    end
end

