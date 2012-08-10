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
        zk_s = [];
        Pk_s = [];
        PPk_s = [];
        logL= [];
    end
    
    methods
        function obj = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL)
            obj.zk_s = zk_s;
            obj.Pk_s = Pk_s;
            obj.PPk_s = PPk_s;
            obj.logL = logL;
        end
        
        function plot_indicator(obj,type,level)
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
            figure
            if strcmp(type,'max')
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