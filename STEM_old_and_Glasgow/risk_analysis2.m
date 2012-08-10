%plot convoluzione e popolazione oltresoglia

clc
load ../Data/krig_mask_scotlandonly2009.mat
load ../Data/popcount_scotland2008.mat
mask(isnotnan(mask))=1;
population=sub;
clear sub
population=population.*mask;

if 1
    load c:\Francesco\STEM3_output\kriging\st_krig_result_NO2_back.mat
    NO2=st_krig_result_NO2.Y_hat;
    NO2_thr=NO2;
    NO2_thr(NO2_thr<105)=0;
    NO2_thr(NO2_thr>=105)=1;
    for t=1:365
        NO2_thr(:,:,t)=NO2_thr(:,:,t).*population;
        NO2_thr_pop(t)=nansum(nansum(NO2_thr(:,:,t)));
    end
    NO2=mean(NO2,3);
    NO2=NO2.*mask;
    data=zeros(nansum(population(:)),1);
    index=1;
    for j=1:size(NO2,2)
        for i=1:size(NO2,1)
            if isnotnan(population(i,j))
                chunk=repmat(NO2(i,j),[population(i,j),1]);
                data(index:index+length(chunk)-1)=chunk;
                index=index+length(chunk);
            end
        end
        j
    end
    [f_NO2,x_NO2]=ksdensity(data,'npoints',1000);
    clear st_krig_result_NO2
    clear NO2_thr
    clear NO2
    
    load c:\Francesco\STEM3_output\kriging\st_krig_result_O3_back.mat
    O3=st_krig_result_O3.Y_hat;
    O3_thr=O3;
    O3_thr(O3_thr<87)=0;
    O3_thr(O3_thr>=87)=1;
    for t=1:365
        O3_thr(:,:,t)=O3_thr(:,:,t).*population;
        O3_thr_pop(t)=nansum(nansum(O3_thr(:,:,t)));
    end
    
    O3=mean(O3,3);
    O3=O3.*mask;
    data=zeros(nansum(population(:)),1);
    index=1;
    for j=1:size(O3,2)
        for i=1:size(O3,1)
            if isnotnan(population(i,j))
                chunk=repmat(O3(i,j),[population(i,j),1]);
                data(index:index+length(chunk)-1)=chunk;
                index=index+length(chunk);
            end
        end
        j
    end
    [f_O3,x_O3]=ksdensity(data,'npoints',1000);
    clear st_krig_result_O3
    clear O3_thr
    clear O3
    
    load c:\Francesco\STEM3_output\kriging\st_krig_result_PM10_back.mat
    PM10=st_krig_result_PM10.Y_hat;
    PM10_thr=PM10;
    PM10_thr(PM10_thr<50)=0;
    PM10_thr(PM10_thr>=50)=1;
    for t=1:365
        PM10_thr(:,:,t)=PM10_thr(:,:,t).*population;
        PM10_thr_pop(t)=nansum(nansum(PM10_thr(:,:,t)));
    end
    
    PM10=mean(PM10,3);
    PM10=PM10.*mask;
    data=zeros(nansum(population(:)),1);
    index=1;
    for j=1:size(PM10,2)
        for i=1:size(PM10,1)
            if isnotnan(population(i,j))
                chunk=repmat(PM10(i,j),[population(i,j),1]);
                data(index:index+length(chunk)-1)=chunk;
                index=index+length(chunk);
            end
        end
        j
    end
    [f_PM10,x_PM10]=ksdensity(data,'npoints',1000);
    clear st_krig_result_PM10
    clear PM10_thr
    clear PM10
    
    figure
    plot(x_NO2,f_NO2,'-','LineWidth',2);
    hold on
    plot(x_O3,f_O3,'--');
    plot(x_PM10,f_PM10,'-')
    figure
    plot((NO2_thr_pop),'s')
    hold on
    plot((O3_thr_pop),'*')
    plot((PM10_thr_pop),'o')
end