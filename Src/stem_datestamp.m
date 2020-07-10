%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
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

classdef stem_datestamp < handle
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    
    %CONSTANTS
    %T     - number of time steps
    
    properties
        date_start=[];      %[integer>0]  (1x1) the date of the first time step. It can be a date in the Matlab format (the output of the datenum function) or a numeric index
        date_end=[];        %[integer>0]  (1x1) the date of the last  time step. It can be a date in the Matlab format (the output of the datenum function) or a numeric index
        T=[];               %[integer>0]  (1x1) the total number of time steps
        step_unit=[]        %[string]     (1x1) unit of time step can be "month", "year", "day", "hour", "minute"
    end
    
    properties (SetAccess = private)
        stamp=[];           %[integer >0] (Tx1) all the date stamps
        irregular=[];       %[integer]    (1x1) 0: time steps are regular, 1: time steps are irregular
    end
    
    methods
        function obj = stem_datestamp(varargin)
            %DESCRIPTION: object constructor
            %
            %INPUT - CASE 1
            %date_start    - [string|datetime]      (1x1) time related to the first time step. It can be a string in the format dd-MM-yyyy HH:mm or a datetime
            %date_end      - [string|datetime]      (1x1) time related to the last  time step. It can be a string in the format dd-MM-yyyy HH:mm or a datetime
            %T             - [integer>0]            (1x1) the total number of time steps
            %<step_unit>   - [string]               (1x1) unit of time step can be "month", "year", "day", "hour", "minute"
            %
            %INPUT - CASE 2
            %dates         - [string]               (Tx1) times related to the T observations. It can be a cell array of string in the format dd-MM-yyyy HH:mm
            %
            %OUTPUT
            %obj           - [stem_datestamp object](1x1)  
            
            if not(nargin==1||nargin==3||nargin==4)
                error('Input arguments must be either 1 or 3 or 4');
            end
            
            if (nargin==1)
                if not(isvector(varargin{1})||iscell(varargin{1}))
                    error('The input argument must be a vector of continuous times or a cell array of dates in the format dd-MM-yyyy HH:mm');
                end
                if (isvector(varargin{1}))
                   if not(isrow(varargin{1}))
                       varargin{1}=varargin{1}';
                   end
                   obj.stamp=varargin{1};
                else
                   obj.stamp=zeros(length(varargin{1}),1);
                   dates=varargin{1};
                   for i=1:length(obj.stamp)
                       obj.stamp(i)=datenum(dates{i},'dd-MM-yyyy HH:mm');
                   end
                end
                if min(diff(obj.stamp))<0
                    error('Times must be ordered ascending');
                end
                obj.date_start=obj.stamp(1);
                obj.date_end=obj.stamp(end);
                obj.T=length(varargin{1});
                obj.irregular=1;
            else
                date_start=varargin{1};
                date_end=varargin{2};

                if isdatetime(date_start)
                    obj.date_start=date_start;
                elseif ischar(date_start)
                    obj.date_start=datetime(date_start,'InputFormat','dd-MM-yyyy HH:mm');
                elseif isstring(date_start)
                    obj.date_start=datetime(date_start,'InputFormat','dd-MM-yyyy HH:mm');
                else
                    error('Please format the date as dd-MM-yyyy HH:mm')
                end

                if isdatetime(date_end)
                    obj.date_end=date_end;
                elseif ischar(date_end)
                    obj.date_end=datetime(date_end,'InputFormat','dd-MM-yyyy HH:mm');
                elseif isstring(date_end)
                    obj.date_end=datetime(date_end,'InputFormat','dd-MM-yyyy HH:mm');
                else
                    error('Please format the date as dd-MM-yyyy HH:mm')
                end

                date_start_num=datenum(obj.date_start);
                date_end_num=datenum(obj.date_end);
                if date_end_num<date_start_num
                    error('date_start cannot be higher than date_end');
                end

                T=varargin{3};
                obj.T=T;
                obj.stamp=obj.date_start:(obj.date_end-obj.date_start)/(obj.T-1):obj.date_end;
                obj.irregular=0;
                if (nargin==4)
                    obj.step_unit=varargin{4};
                    
                    if sum(strcmp(obj.step_unit,["month","year","day","hour","minute"]))==0
                        error('Unit of time step must be "month", "year", "day", "hour", "minute"')
                    end
                    %according to the unit, change the stamp if it is
                    %monthly or yearly.
                    if strcmp(obj.step_unit,'month')
                        obj.stamp=obj.date_start+calmonths(0:T-1);
                    end
                    if strcmp(obj.step_unit,'year')
                        obj.stamp=obj.date_start+calyears(0:T-1);
                    end
                end
                
            end
             
        end
        
        function subset_stamps(obj,indices)
            %DESCRIPTION: removes a subset of the date stamps
            %
            %INPUT
            %obj           - [stem_datestamp object]    (1x1)  stem_datestamp object
            %indices       - [integer>0]                (dTx1) the indices of the temporal steps to keep
            %
            %OUTPUT
            %none: the properties of the object are updated
            
            obj.stamp=obj.stamp(indices);
            obj.date_start=min(obj.stamp);
            obj.date_end=max(obj.stamp);
            obj.T=length(obj.stamp);
        end
        
        function average_stamps(obj,n_steps)
            %DESCRIPTION: averages the date stamps. This method is used when the time_average method of the class stem_data is called
            %
            %INPUT
            %obj           - [stem_datestamp object]    (1x1) stem_datestamp object
            %n_steps       - [integer >0]               (1x1) the number of temporal steps to average            %
            %
            %OUTPUT
            %none: the properties of the object are updated  
            
            indices=0:n_steps:obj.T;
            
            for i=1:length(indices)-1
                %stamp_temp(i)=mean(obj.stamp(indices(i)+1):obj.stamp(indices(i+1)));
                stamp_temp(i)=obj.stamp(indices(i)+1)+...
                    (obj.stamp(indices(i+1))-obj.stamp(indices(i)+1))/2;
            end
            obj.stamp=stamp_temp;
            obj.date_start=min(obj.stamp);
            obj.date_end=max(obj.stamp);
            obj.T=length(obj.stamp);
        end
    end
end