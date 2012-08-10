% clc
% clear all
% 
% load ../Data/st_krig_result_PM10_Glasgow2,mat
% load ../Data/st_loo_residual_Glasgow2.mat
% load ../Data/st_model_Glasgow2.mat
% 
% %recover data transformation details
% log_transformed=st_model.stem_data.stem_varset.log_transformed;
% standardized=st_model.stem_data.stem_varset.standardized;
% means=st_model.stem_data.stem_varset.means;
% stds=st_model.stem_data.stem_varset.stds;
% clear st_model
% %recover loo residual and kriging variances
% resPM10=st_loo_residual.Y(end-60+1:end,:);
% kVarPM10=st_loo_residual.kriging_Var_W_bar_hat(end-60+1:end,:);
% clear st_loo_residual
% %studentized loo residuals
% resPM10_studentized=resPM10./sqrt(kVarPM10);
% 
% %set the threshold
% thrPM10=50;
% original_thr=thrPM10;
% %transform the threshold if necessary
% if log_transformed
%     thrPM10=log(thrPM10);
% end
% if standardized
%     thrPM10=(thrPM10-means(3))/stds(3);
% end
% clc
% %evaluate the cdf of the studentized residuals
% [fPM10,xPM10]=ksdensity(resPM10_studentized(:),'npoints',1000,'function','cdf');
% %plot(xPM10,fPM10);
% 
% 
% for t=1:365
%     t
%     files=dir(['../Data/multi_krig/t',num2str(t),'/*.mat']);
%     thr_prob=zeros(size(st_krig_result_PM10.Y_hat,1),size(st_krig_result_PM10.Y_hat,2),length(files));
%     for b=1:1%length(files)
%         load(['../Data/multi_krig/t',num2str(t),'/',files(b).name]);
%         for j=1:size(st_krig_result_PM10.Y_hat,2)
%             for i=1:size(st_krig_result_PM10.Y_hat,1)
%                 if isnotnan(st_krig_result_PM10.Y_hat(i,j,t))
%                     x=(xPM10*sqrt(st_krig_result_PM10.Var_Y_hat(i,j,t))) + boot_krig_fixtime.Y_hat(i,j);
%                     if thrPM10>max(x)
%                         thr_prob(i,j,b)=0;
%                     else
%                         if thrPM10<min(x)
%                             thr_prob(i,j,b)=1;
%                         else
%                             L=thrPM10<x;
%                             index=find(L,1);
%                             thr_prob(i,j,b)=1-fPM10(index);
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     save(['../Data/multi_risk/risk_t',num2str(t),'_level',num2str(original_thr)],'thr_prob','-v7.3');
% end

clc
clear all
load ../Data/popcount_scotland2008.mat
load ../Data/krig_mask_scotlandonly2009.mat
mask(isnotnan(mask))=1;
sub=sub.*mask;
files=dir(['../Data/multi_risk/*.mat']);
index=zeros(365,1);
for t=1:length(files)
    t
    load(['../Data/multi_risk/',files(t).name])
    name=files(t).name;
    code=name(7:end-12);
    temp=thr_prob(:,:,1);
    temp=temp.*sub;
    index(str2num(code))=nansum(temp(:));
end





