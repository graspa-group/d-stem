%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_kalmansmoother_result < handle
    properties
        zk_s  = [];     %[double]     (pxT+1)    the smoothed state
        Pk_s  = [];     %[double]     (pxpxT+1)  variance-covariance matrix of the smoothed state
        PPk_s = [];     %[double]     (pxpxT+1)  lag-one variance-covariance matrix of the smoothed state
        logL  = [];     %[double]     (1x1)      observed-data log-likelihood
    end
    
    methods
        function obj = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL)
            %DESCRIPTION: constructor of the class stem_kalmansmoother_result
            %
            %INPUT 
            %See the class properties
            %
            %OUTPUT
            %obj             - [stem_kalmansmoother_result object]   (1x1) stem_kalmansmoother_result object
            
            obj.zk_s = zk_s;
            obj.Pk_s = Pk_s;
            obj.PPk_s = PPk_s;
            obj.logL = logL;
        end
        
        function plot(obj,level,flag_max)
            %DESCRIPTION: plot zk_s and confidence bands based on Pk_s
            %
            %INPUT 
            %obj            - [stem_kalmansmoother_result object]   (1x1) stem_kalmansmoother_result object
            %<level>        - [double >0 and <1] (default: 0.05)    (1x1) the confidence level
            %<flag_max>     - [boolean]                             (1x1) 1: also max(zk_s) is plotted; 0: no additional plot
            %
            %OUTPUT
            %none: zk_s is plotted
            
            if nargin<2
                level=0.05;
            end
            if nargin<3
                flag_max=0;
            end
              
            marker_list={'x','o','^','s'};
            l=size(obj.zk_s,1);
            if l<=3
                rows=l;
                cols=1;
            else
                l=size(obj.zk_s,1)^0.5;
                if round(l)^2==size(obj.zk_s,1)
                    rows=l;
                    cols=l;
                else
                    rows=ceil(l);
                    cols=round(l);
                end
            end
            
            figure
            for i=1:size(obj.zk_s,1)
                subplot(rows,cols,i);
                plot(obj.zk_s(i,2:end)');
                hold on
                plot(obj.zk_s(i,2:end)'+2*squeeze(sqrt(obj.Pk_s(i,i,2:end))),':');
                plot(obj.zk_s(i,2:end)'-2*squeeze(sqrt(obj.Pk_s(i,i,2:end))),':');
            end
            if (flag_max)
                figure
                for t=2:length(obj.zk_s)
                    t
                    obj.Pk_s(1,2,t)=obj.Pk_s(2,1,t);
                    obj.Pk_s(1,3,t)=obj.Pk_s(3,1,t);
                    obj.Pk_s(2,3,t)=obj.Pk_s(3,2,t);
                    ma=max(obj.zk_s(:,t));
                    step=ma/1000;
                    p=0.5;
                    while (p<1-level/2)
                        ma=ma+step;
                        p=mvncdf(repmat(ma,length(obj.zk_s(:,t)),1),obj.zk_s(:,t),obj.Pk_s(:,:,t));
                    end
                    up(t-1)=ma;
                    ma=max(obj.zk_s(:,t));
                    p=0.5;
                    while (p>level/2)
                        ma=ma-step;
                        p=mvncdf(repmat(ma,length(obj.zk_s(:,t)),1),obj.zk_s(:,t),obj.Pk_s(:,:,t));
                    end
                    lo(t-1)=ma;
                end
                [val,idx]=max(obj.zk_s(:,2:end));
                hold on
                for t=1:length(val)
                    plot(t,val(t),'Marker',marker_list{idx(t)},'MarkerSize',5);
                    plot([t,t],[lo(t),up(t)]);
                end
            end
        end
    end
end