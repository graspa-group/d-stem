clc
clear all
load ../Output/Paper/no2_hourly.mat
load ../Output/Paper/o3_hourly.mat

% counter=1;
% for i=1:67
%     i
%     for t=1:365
%         data=[]; 
%         for k=1:17
%             temp=no2(i,(t-1)*24+k:(t-1)*24+k+7);
%             data=[data nanmean(temp)];
%         end
%         running_max(counter)=max(data);
%         daily_mean(counter)=nanmean(no2(i,(t-1)*24+1:t*24));
%         counter=counter+1;
%     end
% end
% plot(daily_mean,running_max,'.')

counter=1;
for i=1:11
    for t=1:365
        data=[]; 
        for k=1:17
            temp=o3(i,(t-1)*24+k:(t-1)*24+k+7);
            data=[data nanmean(temp)];
        end
        running_max(counter)=max(data);
        daily_mean(counter)=nanmean(o3(i,(t-1)*24+1:t*24));
        counter=counter+1;
    end
end
subplot(1,2,2);
L1=isnotnan(daily_mean);
L2=isnotnan(running_max);
L=L1&L2;
daily_mean=daily_mean(L);
running_max=running_max(L);
plot(daily_mean,running_max,'.')
hold on
a=inv(daily_mean*daily_mean')*daily_mean*running_max'
x=min(daily_mean):0.1:max(daily_mean);
y=x*a;
plot(x,y,'r-');
axis equal

clear running_max
clear daily_mean

counter=1;
for i=1:67
    for t=1:365
        running_max(counter)=max(no2(i,(t-1)*24+1:t*24));
        daily_mean(counter)=nanmean(no2(i,(t-1)*24+1:t*24));
        counter=counter+1;
    end
end
subplot(1,2,1);
L1=isnotnan(daily_mean);
L2=isnotnan(running_max);
L=L1&L2;
daily_mean=daily_mean(L);
running_max=running_max(L);
plot(daily_mean,running_max,'.')
hold on
a=inv(daily_mean*daily_mean')*daily_mean*running_max'
x=min(daily_mean):0.1:max(daily_mean);
y=x*a;
plot(x,y,'r-');
axis equal
