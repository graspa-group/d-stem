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

classdef stem_krig_options
    properties
        block_size=0;               %[integer>0]        (1x1)   the size of the kriging blocks. If set to 0, only one block is used.
        nn_size=0;                  %[integer>0]        (1x1)   the size of the nearest neighbour set for each kriging location. If set to 0, the entire dataset is used
        back_transform=1;           %[boolean]          (1x1)   1: kriging output are back-transformed to the orginal unit; 0: kriging output are NOT back transformed 
        no_varcov=1;                %[boolean]          (1x1)   1: the variance of the kriging output is NOT computed; 0: the variance is computed
        crossval=0;                 %[boolean]          (1x1)   1: the kriging is done for cross-validation; 0: standard kriging. This flag must be 0 when you are doing kriging
        type='y';                   %[string]           (1x1)   'y': the kriging is on the y variable(s), 'y-xbeta': the kriging is on y-X_beta*beta', this option is only valid with f-HDGM models
        workers=1;                  %[integer>0]        (1x1)   the number of matlab workers used for block kriging
    end
    
    methods
        function obj = stem_krig_options()
            %DESCRIPTION: object constructor
            %
            %INPUT
            %no inputs required
            %
            %OUTPUT
            %obj       - [stem_krig_options] (1x1)
        end
        
        %Class set methods
        function obj = set.block_size(obj,block_size)
            if (block_size<0)
                error('block_size must be >=0');
            end
            obj.block_size=block_size;
        end
        
        function obj = set.nn_size(obj,nn_size)
            if (nn_size<0)
                error('nn_size must be >=0');
            end
            obj.nn_size=nn_size;
        end
        
        function obj = set.back_transform(obj,back_transform)
            if not(back_transform==0||back_transform==1)
                error('back_trasform must be either 0 or 1');
            end
            obj.back_transform=back_transform;
        end
        
        function obj = set.no_varcov(obj,no_varcov)
            if not(no_varcov==0||no_varcov==1)
                error('no_varcov must be either 0 or 1');
            end
            obj.no_varcov=no_varcov;
        end
        
        function obj = set.crossval(obj,crossval)
            if not(crossval==0||crossval==1)
                error('crossval must be either 0 or 1');
            end
            obj.crossval=crossval;
        end
        
        function obj = set.type(obj,type)
            if not(strcmp(type,'y')||strcmp(type,'y-xbeta'))
                error('type must be either ''y'' or ''y-xbeta''');
            end
            obj.type=type;
        end
        
        function obj = set.workers(obj,workers)
            if workers<=0
                error('workers must be >0');
            end
            obj.workers=workers;
        end

    end
end