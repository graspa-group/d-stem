clc
clear all
load ../Data/population.mat
pop2(pop2==0)=NaN;
load ../Data/st_krig_result__20120603_135334
T=size(st_krig_result.y_hat,3);
for t=1:T
    E=zeros(nansum(pop2(:)),1);
    idx=0;
    for j=1:size(st_krig_result.y_hat,1)
        j
        for i=1:size(st_krig_result.y_hat,2)
            if not(isnan(pop2(j,i)))
                v=repmat(st_krig_result.y_hat(j,i,t),pop2(j,i),1);
                l=length(v);
                E(idx+1:idx+l,1)=v;
                idx=idx+l;
            end
        end
    end
    [y,x]=hist(E,50);
    Y(:,t)=y;
    X(:,t)=x;
end
for i=1:35
    color = i/size(Y,2);
    c=[0 color color];
    scs = csaps(0:5:300,Y(:,i));
    y=fnval(scs,0:0.1:300);
    plot(0:0.1:300,y,'Color',c,'LineWidth',2);
    hold on
end

