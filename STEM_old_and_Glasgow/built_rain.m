clc
clear all
%carica la mappa di destinazione

if 0
    lat_start=43;
    lat_end=47;
    lon_start=6;
    lon_end=14;
    step=0.25;
    lat=[lat_start:step:lat_end]';
    lon=[lon_start:step:lon_end]';
    rain=zeros(length(lat),length(lon),365*8);
    %leggo i nomi dei file nella cartella dei dati
    cartella='C:\Documents and Settings\Francesco\Desktop\Rain\*.hdf';
    files=dir(cartella);
    day_index=1;
    counter=0;
    for i=1:length(files)
        %estrae separatamente ora,giorno,mese e anno
        day=files(i).name(10:11);
        month=files(i).name(8:9);
        year=files(i).name(6:7);
        hour=files(i).name(13:14);
        h=str2num(hour)/3;
        date=datestr([month '/' day '/' year]);
        %converte la data in numero
        d=datenum(date);
        if i==1
            first_day=d;
        end
        fname=[cartella(1:end-5) files(i).name];
        data = hdfread(fname,'precipitation');
        data=rot90(data);
        lat_index=round((50-lat_end)/0.25):1:round((50-lat_start)/0.25);
        lon_index=round((180+lon_start)/0.25):1:round((180+lon_end)/0.25);
        data=double(data(lat_index,lon_index));
        rain(:,:,(d-first_day)*8+1+h)=data;
        i
    end
    save rain_hourly_2006 rain
end

if 0
    clc
    load rain_hourly_2006.mat
    s=size(rain);
    rain_daily=zeros(s(1),s(2),s(3)/8);
    for i=1:365
     rain_daily(:,:,i)=sum(rain(:,:,(i-1)*8+1:i*8),3)*3;    
    end
    save rain_daily_2006 rain_daily
    save rain_lat lat
    save rain_lon lon
    figure
end

if 1
    load rain_daily_2006
    for i=1:365
        temp=rain_daily(:,:,i);
        temp(1,1)=max(rain_daily(:));
        temp(1,2)=0;
        imagesc(temp);
        title(['Day: ' num2str(i)]);
        colorbar
        drawnow
        pause(0.3);
    end
end



