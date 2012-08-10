clc
load ../Data/popcount_scotland2008.mat
load ../Data/krig_mask_scotlandonly2009.mat
mask(isnotnan(mask))=1;
sub=sub.*mask;
m=[31 28 31 30 31 30 31 31 30 31 30 31];
block=[0 cumsum(m)];
for i=1:length(m)
    y_NO2=mean(st_krig_result_NO2.Y_hat(:,:,blocks(i)+1:blocks(i+1)),3);
    NO2_indicator(i)=nansum(nansum(y_NO2.*sub))/nansum(sub(:));
end
for i=1:length(m)
    y_O3=mean(st_krig_result_O3.Y_hat(:,:,blocks(i)+1:blocks(i+1)),3);
    O3_indicator(i)=nansum(nansum(y_O3.*sub))/nansum(sub(:));
end
for i=1:length(m)
    y_pm10=mean(st_krig_result_PM10.Y_hat(:,:,blocks(i)+1:blocks(i+1)),3);
    pm10_indicator(i)=nansum(nansum(y_pm10.*sub))/nansum(sub(:));
end

subplot(3,1,1);
plot(NO2_indicator,'*-')
hold on
plot(ones(12,1)*mean(NO2_indicator),'--');
subplot(3,1,2);
plot(O3_indicator,'*-')
hold on
plot(ones(12,1)*mean(O3_indicator),'--');
subplot(3,1,3);
plot(pm10_indicator,'*-')
hold on
plot(ones(12,1)*mean(pm10_indicator),'--');
