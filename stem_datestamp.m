%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_datestamp < handle
    properties
        date_start=[];
        date_end=[];
        T=[];
    end
    
    properties (SetAccess = private)
        stamp=[];
    end
    
    methods
        function obj = stem_datestamp(date_start,date_end,T)
            if nargin<3
                error('Not enough input arguments');
            end
            if isnumeric(date_start)
                obj.date_start=date_start;
            else
                obj.date_start=datenum(date_start,'dd-mm-yyyy');
            end
            if isnumeric(date_end)
                obj.date_end=date_end;
            else
                obj.date_end=datenum(date_end,'dd-mm-yyyy');
            end
            if obj.date_end<obj.date_start
                error('date_start cannot be higher than date_end');
            end
            obj.T=T;
            obj.stamp=obj.date_start:(obj.date_end-obj.date_start+1)/obj.T:obj.date_end;
        end
        
        function subset_stamps(obj,indices)
            obj.stamp=obj.stamp(indices);
            obj.date_start=min(obj.stamp);
            obj.date_end=max(obj.stamp);
            obj.T=length(obj.stamp);
        end
    end
    
    methods (Static)
        function res = compare(stem_datestamp1,stem_datestamp2)
            if sum(stem_datestamp1.stamp-stem_datestamp2.stamp)==0
                res=1;
            else
                res=0;
            end
        end
    end
end