function [raw_data] = build_raw_data(variables,covariates)
disp('Building raw data...')
if nargin<1
    flag_const=0;
end
T=365;

variables_counter=1;
flag_omi=1;

%%%%%%%%%%%%%%%%%%%%%%%% AOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('aot',variables))
    
    if flag_omi
        load ../Data/AOT_2006/aot_omi.mat
        aot_omi_dest=zeros(54,21,T);
        for t=1:T
            aot_omi_dest(:,:,t)=fliplr(aot_omi(:,:,t)');
        end
        
        aot_omi_lat=fliplr(aot_omi_lat');
        aot_omi_long=aot_omi_long';
        aot_lat=reshape(aot_omi_lat,54*21,1);
        aot_long=reshape(aot_omi_long,54*21,1);
        aot_data=reshape(aot_omi_dest,54*21,T);
        aot_coordinates=[aot_lat aot_long];
    else
        load ../Data/AOT_2006/aot_2006.mat
        load ../Data/AOT_2006/aot_lat.mat
        load ../Data/AOT_2006/aot_lon.mat
        
        lat_dest=fliplr(lat_dest);
        for i=1:T
            temp=aot_data(:,:,i);
            temp=fliplr(temp);
            aot_data_flip(:,:,i)=temp;
        end
        aot_data=aot_data_flip;
        
        aot_lat=reshape(lat_dest,54*21,1);
        aot_long=reshape(long_dest,54*21,1);
        aot_coordinates=[aot_lat aot_long];
        aot_data=reshape(aot_data,54*21,T);
    end
    
    raw_data.Y{variables_counter}.data=aot_data;
    raw_data.Y{variables_counter}.grid=aot_coordinates;
    raw_data.Y{variables_counter}.grid_size=[54 21];
    raw_data.Y{variables_counter}.name='aot';
    raw_data.Y{variables_counter}.grid_type='regular';
    
    variables_counter=variables_counter+1;
end

%%%%%%%%%%%%%%%%%%%%%%% PM10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('pm10',variables))
    pm_data_lombardia = csvread('../Data/PM10/PM10_lombardia2006.csv')';
    pm_data_piemonte = csvread('../Data/PM10/PM10_piemonte2006.csv')';
    pm_data_veneto = csvread('../Data/PM10/PM10_veneto2006.csv')';
    pm_coordinates_lombardia = csvread('../Data/PM10/anagrafica_Lombardia.csv');
    pm_coordinates_piemonte = csvread('../Data/PM10/anagrafica_Piemonte.csv');
    pm_coordinates_veneto = csvread('../Data/PM10/anagrafica_Veneto.csv');
    [anagrafica_pianura,txt,raw] = xlsread('../Data/PM10/anagrafica_PianuraPadana.xls');
    
    pm_code_lombardia=pm_coordinates_lombardia(:,1);
    pm_coordinates_lombardia=pm_coordinates_lombardia(:,2:end);
    pm_code_piemonte=pm_coordinates_piemonte(:,1);
    pm_coordinates_piemonte=pm_coordinates_piemonte(:,2:end);
    pm_code_veneto=pm_coordinates_veneto(:,1);
    pm_coordinates_veneto=pm_coordinates_veneto(:,2:end);
    pm_data_lombardia=pm_data_lombardia(:,2:end);
    pm_data_piemonte=pm_data_piemonte(:,2:end);
    pm_data_veneto=pm_data_veneto(:,2:end);
    
    pm_data=[pm_data_lombardia;pm_data_piemonte;pm_data_veneto];
    pm_code=[pm_code_lombardia;pm_code_piemonte;pm_code_veneto];
    pm_coordinates=[pm_coordinates_lombardia;pm_coordinates_piemonte;pm_coordinates_veneto];
    pm_coordinates = circshift(pm_coordinates,[0 1]);
    
    for i=1:length(pm_code)
        ref=pm_code(i);
        for j=1:length(anagrafica_pianura)
            if anagrafica_pianura(j,1)==ref
                pm_name{i}=txt(j,4);
            end
        end
    end
    
    %filtro le centraline PM rispetto ad una soglia di missing
    max_missing_rate=0.3;
    miss_rate=sum(isnan(pm_data)')/T;
    pm_data=pm_data';
    pm_data=pm_data(:,miss_rate<max_missing_rate)';
    pm_coordinates=pm_coordinates';
    pm_coordinates=pm_coordinates(:,miss_rate<max_missing_rate)';
    pm_name=pm_name(miss_rate<max_missing_rate);
    
    raw_data.Y{variables_counter}.data=pm_data;
    raw_data.Y{variables_counter}.grid=pm_coordinates;
    raw_data.Y{variables_counter}.name=pm_name;
    raw_data.Y{variables_counter}.grid_size=[];
    raw_data.Y{variables_counter}.name='pm';
    raw_data.Y{variables_counter}.grid_type='sparse';
end

covariates_counter=1;

%%%%%%%%%%%%%%%%%%%%%%%% BLH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('blh',covariates))
    blh1=load('-ascii','../Data/BLH_2006/blh_01_06.ascii');
    blh2=load('-ascii','../Data/BLH_2006/blh_02_06.ascii');
    blh3=load('-ascii','../Data/BLH_2006/blh_03_06.ascii');
    blh4=load('-ascii','../Data/BLH_2006/blh_04_06.ascii');
    blh5=load('-ascii','../Data/BLH_2006/blh_05_06.ascii');
    blh6=load('-ascii','../Data/BLH_2006/blh_06_06.ascii');
    blh7=load('-ascii','../Data/BLH_2006/blh_07_06.ascii');
    blh8=load('-ascii','../Data/BLH_2006/blh_08_06.ascii');
    blh9=load('-ascii','../Data/BLH_2006/blh_09_06.ascii');
    blh10=load('-ascii','../Data/BLH_2006/blh_10_06.ascii');
    blh11=load('-ascii','../Data/BLH_2006/blh_11_06.ascii');
    blh12=load('-ascii','../Data/BLH_2006/blh_12_06.ascii');
    blh = [blh1;blh2;blh3;blh4;blh5;blh6;blh7;blh8;blh9;blh10;blh11;blh12];
    clear blh1 blh2 blh3 blh4 blh5 blh6 blh7 blh8 blh9 blh10 blh11 blh12;
    blh_data = reshape(blh(:,3),1681,T);
    blh_coordinates = blh(1:1681,1:2);
    
    raw_data.X{covariates_counter}.data=blh_data;
    raw_data.X{covariates_counter}.grid=blh_coordinates;
    raw_data.X{covariates_counter}.grid_size=[41 41];
    raw_data.X{covariates_counter}.name='blh';
    
    covariates_counter=covariates_counter+1;
end


%%%%%%%%%%%%%%%%%%%%%%%% Elevation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('elevation',covariates))
    elevation=load('-ascii','../Data/Elevation/elevation.ascii');
    long_X=reshape(elevation(:,1),1239,426);
    lat_X=reshape(elevation(:,2),1239,426);
    data_X=reshape(elevation(:,3),1239,426);
    data_X(data_X==-9999)=0;
    
    long_X=reshape(long_X(3:3:1239,3:3:426),413*142,1);
    lat_X=reshape(lat_X(3:3:1239,3:3:426),413*142,1);
    data_X=reshape(data_X(3:3:1239,3:3:426),413*142,1);
    
    elevation_data=repmat(data_X,1,T);
    elevation_coordinates=[lat_X long_X];
    
    raw_data.X{covariates_counter}.data=elevation_data;
    raw_data.X{covariates_counter}.grid=elevation_coordinates;
    raw_data.X{covariates_counter}.grid_size=[413 142];
    raw_data.X{covariates_counter}.name='elevation';
    
    covariates_counter=covariates_counter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%% LAND USE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('urban',covariates))
    load ../Data/Land/land_use_2006.mat
    load ../Data/Land/land_lat.mat
    load ../Data/Land/land_lon.mat
    
    lat=lat';
    lon=lon';
    
    lat=repmat(fliplr(lat'),161,1);
    lon=repmat(lon,1,81);
    
    lat=reshape(lat,161*81,1);
    lon=reshape(lon,161*81,1);
    
    land_coordinates=[lat lon];
    
    for i=1:365
        temp=land_use(:,:,i);
        temp=rot90(temp);
        temp=flipud(temp);
        land_use_r(:,i)=reshape(temp,161*81,1);
    end
    
    raw_data.X{covariates_counter}.data=land_use_r;
    raw_data.X{covariates_counter}.grid=land_coordinates;
    raw_data.X{covariates_counter}.grid_size=[161 81];
    raw_data.X{covariates_counter}.name='urban_density';
    
    covariates_counter=covariates_counter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%% RAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(strcmp('rain',covariates))
    load ../Data/Rain/rain_daily_2006.mat
    load ../Data/Rain/rain_lat.mat
    load ../Data/Rain/rain_lon.mat
    
    lat=repmat(fliplr(lat'),33,1);
    lon=repmat(lon,1,17);
    
    lat=reshape(lat,33*17,1);
    lon=reshape(lon,33*17,1);
    
    rain_coordinates=[lat lon];
    
    for i=1:365
        temp=rain_daily(:,:,i);
        temp=rot90(temp);
        temp=flipud(temp);
        rain_r(:,i)=reshape(temp,33*17,1);
    end
    
    raw_data.X{covariates_counter}.data=rain_r;
    raw_data.X{covariates_counter}.grid=rain_coordinates;
    raw_data.X{covariates_counter}.grid_size=[33 17];
    raw_data.X{covariates_counter}.name='rain';
    
    covariates_counter=covariates_counter+1;
end

