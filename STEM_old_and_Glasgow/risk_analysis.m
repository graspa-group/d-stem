clc
clear all

if 0
    load ../Output/Risk/st_krig_result_PM10
    load ../Output/Risk/st_loo_residual.mat
    load ../Output/Risk/st_model_afterlastloo.mat
    %recover data transformation details
    log_transformed=st_model.stem_data.stem_varset.log_transformed;
    standardized=st_model.stem_data.stem_varset.standardized;
    means=st_model.stem_data.stem_varset.means;
    stds=st_model.stem_data.stem_varset.stds;
    %recover loo residual and kriging variances
    resPM10=st_loo_residual.Y(end-60+1:end,:);
    kVarPM10=st_loo_residual.kriging_Var_W_bar_hat(end-60+1:end,:);
    %studentized loo residuals
    resPM10_studentized=resPM10./kVarPM10;
    
    %set the threshold
    thrPM10=50;
    %transform the threshold if necessary
    if log_transformed
        thrPM10=log(thrPM10);
    end
    if standardized
        thrPM10=(thrPM10-means(3))/stds(3);
    end
    clc
    %evaluate the icdf of the studentized residuals
    [fPM10,xPM10]=ksdensity(resPM10_studentized(:),'npoints',1000,'function','icdf');
    %plot(xPM10,fPM10);
    thr_prob=zeros(size(st_krig_result_PM10.Y_hat));
    %generate random from the distribution (here to save time or at * for a real bootstrap)
    data_studentized = rnd_from_icdf_ksdensity(fPM10,xPM10,500);
    for t=1:size(st_krig_result_PM10.Y_hat,3)
        disp(num2str(t));
        for j=1:size(st_krig_result_PM10.Y_hat,2)
            for i=1:size(st_krig_result_PM10.Y_hat,1)
                if isnotnan(st_krig_result_PM10.Y_hat(i,j,t))
                    %data_studentized = rnd_from_icdf_ksdensity(fPM10,xPM10,500); (*)
                    
                    %back transform the residual
                    data = (data_studentized*st_krig_result_PM10.Var_Y_hat(i,j,t)) + st_krig_result_PM10.Y_hat(i,j,t);
                    
                    [f,x]=ecdf(data); % or but much slower [f,x]=ksdensity(data,'npoints',200,'function','cdf');
                    if thrPM10>max(x)
                        thr_prob(i,j,t)=0;
                    else
                        if thrPM10<min(x)
                            thr_prob(i,j,t)=0;
                        else
                            L=thrPM10<x;
                            index=find(L,1);
                            thr_prob(i,j,t)=1-f(index);
                        end
                    end
                end
            end
        end
    end
    
    save thr_prob_PM10_higher50 thr_prob -v7.3
end

if 0
    load ../Output/Risk/st_krig_result_NO2
    load ../Output/Risk/st_loo_residual.mat
    load ../Output/Risk/st_model_afterlastloo.mat
    %recover data transformation details
    log_transformed=st_model.stem_data.stem_varset.log_transformed;
    standardized=st_model.stem_data.stem_varset.standardized;
    means=st_model.stem_data.stem_varset.means;
    stds=st_model.stem_data.stem_varset.stds;
    %recover loo residual and kriging variances
    resNO2=st_loo_residual.Y(1:66,:);
    kVarNO2=st_loo_residual.kriging_Var_W_bar_hat(1:66,:);
    %studentized loo residuals
    resNO2_studentized=resNO2./kVarNO2;
    
    %set the threshold
    thrNO2=40;
    %transform the threshold if necessary
    if log_transformed
        thrNO2=log(thrNO2);
    end
    if standardized
        thrNO2=(thrNO2-means(1))/stds(1);
    end
    clc
    %evaluate the icdf of the studentized residuals
    [fNO2,xNO2]=ksdensity(resNO2_studentized(:),'npoints',1000,'function','icdf');
    %plot(xNO2,fNO2);
    thr_prob=zeros(size(st_krig_result_NO2.Y_hat));
    %generate random from the distribution (here to save time or at * for a real bootstrap)
    data_studentized = rnd_from_icdf_ksdensity(fNO2,xNO2,500);
    for t=1:size(st_krig_result_NO2.Y_hat,3)
        disp(num2str(t));
        for j=1:size(st_krig_result_NO2.Y_hat,2)
            for i=1:size(st_krig_result_NO2.Y_hat,1)
                if isnotnan(st_krig_result_NO2.Y_hat(i,j,t))
                    %data_studentized = rnd_from_icdf_ksdensity(fNO2,xNO2,500); (*)
                    
                    %back transform the residual
                    data = (data_studentized*st_krig_result_NO2.Var_Y_hat(i,j,t)) + st_krig_result_NO2.Y_hat(i,j,t);
                    
                    [f,x]=ecdf(data); % or but much slower [f,x]=ksdensity(data,'npoints',200,'function','cdf');
                    if thrNO2>max(x)
                        thr_prob(i,j,t)=0;
                    else
                        if thrNO2<min(x)
                            thr_prob(i,j,t)=0;
                        else
                            L=thrNO2<x;
                            index=find(L,1);
                            thr_prob(i,j,t)=1-f(index);
                        end
                    end
                end
            end
        end
    end
    
    save thr_prob_NO2_higher50 thr_prob -v7.3
end

if 1
    load ../Output/Risk/st_krig_result_O3
    load ../Output/Risk/st_loo_residual.mat
    load ../Output/Risk/st_model_afterlastloo.mat
    %recover data transformation details
    log_transformed=st_model.stem_data.stem_varset.log_transformed;
    standardized=st_model.stem_data.stem_varset.standardized;
    means=st_model.stem_data.stem_varset.means;
    stds=st_model.stem_data.stem_varset.stds;
    %recover loo residual and kriging variances
    resO3=st_loo_residual.Y(67:76,:);
    kVarO3=st_loo_residual.kriging_Var_W_bar_hat(67:76,:);
    %studentized loo residuals
    resO3_studentized=resO3./kVarO3;
    
    %set the threshold
    thrO3=100;
    %transform the threshold if necessary
    if log_transformed
        thrO3=log(thrO3);
    end
    if standardized
        thrO3=(thrO3-means(2))/stds(2);
    end
    clc
    %evaluate the icdf of the studentized residuals
    [fO3,xO3]=ksdensity(resO3_studentized(:),'npoints',1000,'function','icdf');
    %plot(xO3,fO3);
    thr_prob=zeros(size(st_krig_result_O3.Y_hat));
    %generate random from the distribution (here to save time or at * for a real bootstrap)
    data_studentized = rnd_from_icdf_ksdensity(fO3,xO3,500);
    for t=1:size(st_krig_result_O3.Y_hat,3)
        disp(num2str(t));
        for j=1:size(st_krig_result_O3.Y_hat,2)
            for i=1:size(st_krig_result_O3.Y_hat,1)
                if isnotnan(st_krig_result_O3.Y_hat(i,j,t))
                    %data_studentized = rnd_from_icdf_ksdensity(fO3,xO3,500); (*)
                    
                    %back transform the residual
                    data = (data_studentized*st_krig_result_O3.Var_Y_hat(i,j,t)) + st_krig_result_O3.Y_hat(i,j,t);
                    
                    [f,x]=ecdf(data); % or but much slower [f,x]=ksdensity(data,'npoints',200,'function','cdf');
                    if thrO3>max(x)
                        thr_prob(i,j,t)=0;
                    else
                        if thrO3<min(x)
                            thr_prob(i,j,t)=0;
                        else
                            L=thrO3<x;
                            index=find(L,1);
                            thr_prob(i,j,t)=1-f(index);
                        end
                    end
                end
            end
        end
    end
    
    save thr_prob_O3_higher50 thr_prob -v7.3
end
