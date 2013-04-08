%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-STEM - Distributed Space Time Expecation Maximization      %
%                                                              %
% Author: Francesco Finazzi                                    %
% E-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo - Dept. of Engineering    %
% Author website: http://www.unibg.it/pers/?francesco.finazzi  %
% Code website: https://code.google.com/p/d-stem/              %
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

