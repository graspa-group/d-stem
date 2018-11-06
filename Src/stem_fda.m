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

classdef stem_fda
   
    properties (SetAccess=private)
        spline_order=[];    %[integer >0]   (1x1) the order of the spline
        spline_range=[];    %[double]       (2x1) the range of the spline
        spline_knots=[];    %[double]       (sx1) the spline knots
        spline_basis=[];    %[basis obj]    (1x1) the spline basis 
    end
    
    methods
        function obj = stem_fda(spline_order,spline_range,spline_knots)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %spline_order=[]    -   [integer >0]        (1x1) the order of the spline
            %spline_range=[]    -   [double]            (2x1) the range of the spline
            %spline_knots=[]    -   [double]            (sx1) the spline knots
            %
            %OUTPUT
            %obj                -   [stem_fda object]   (1x1)
            
            if nargin<3
                error('All the input arguments must be provided');
            end
            obj.spline_order=spline_order;
            obj.spline_range=spline_range;
            obj.spline_knots=spline_knots;
            
            if max(spline_knots)>spline_range(2)
                error('The highest element of spline_knots cannot be higher than spline_range(2)');
            end
            if min(spline_knots)<spline_range(1)
                error('The lowest element of spline_knots cannot be lower than spline_range(1)');
            end
            
            norder=obj.spline_order+1;
            nbasis=length(obj.spline_knots)+norder-2;
            obj.spline_basis=create_bspline_basis(obj.spline_range, nbasis, norder, obj.spline_knots);
        end
             
        %Class set methods
        function obj = set.spline_order(obj,spline_order)
            if spline_order<1
                error('spline_order must be > 0');
            end
            obj.spline_order=spline_order;
        end
        
        function obj = set.spline_range(obj,spline_range)
            if not(length(spline_range)==2)
                error('spline_range must be a 2x1 vector');
            end
            if spline_range(2)<=spline_range(1)
                error('The upper bound of spline_range must be higher than the lower bound');
            end
            obj.spline_range=spline_range;
        end
        
        function obj = set.spline_knots(obj,spline_knots)
            if length(spline_knots)<2
                error('spline_knots must be at least 2x1');
            end
            if any(diff(spline_knots)<=0)
                error('spline_knots must be a vector of sorted and non equal elements');
            end
            obj.spline_knots=spline_knots;
        end
        
    end
end