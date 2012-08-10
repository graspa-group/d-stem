clc
clear all
%carica la mappa di destinazione
lat_start=44;
lat_end=47;
long_start=6;
long_end=14;
step=0.15;
lat=[lat_start:step:lat_end];
long=[long_start:step:long_end];
[lat_dest,long_dest] = meshgrid(lat,long);
lat_dest=lat_dest';
long_dest=long_dest';
if 1
    %leggo i nomi dei file nella cartella dei dati
    cartella='E:\Ricerca\Stem2\AOT_Raw\2006\';
    files=dir(cartella);
    aot_dest=zeros(size(lat_dest,1),size(lat_dest,2),365,3)*NaN;
    aot_flag=zeros(365,3);
    last_day=0;
    day_index=1;
    counter=0;
    for i=3:length(files)
        day=str2num(files(i).name(15:17));
        if day==last_day
            day_index=day_index+1;
        else
            day_index=1;
            counter=counter+1;
        end
        time=files(i).name(19:22);
        fname=[cartella files(i).name];
        lat = hdfread(fname, 'mod04','Fields','Latitude');
        long = hdfread(fname, 'mod04','Fields','Longitude');
        aot = hdfread(fname, 'mod04','Fields','Optical_Depth_Land_And_Ocean');
        
        %regions = shaperead('..\maps\italy\Ita_Reg.shp');
        %regions(1,1).Geometry='Line';
        aot=double(aot);
        lat=double(lat);
        long=double(long);
        aot(aot==-9999)=NaN;
        %     figure;
        %     mapshow(long,lat,aot,'DisplayType','surface');
        %     hold on
        %     mapshow(regions);
        if size(aot,1)*size(aot,2)>4
            aot_dest(:,:,counter,day_index)=griddata(lat,long,aot,lat_dest,long_dest);
            aot_flag(counter,day_index)=1;
        end
        last_day=day;
        %     figure;
        %     mapshow(long_dest,lat_dest,aot_dest,'DisplayType','surface');
        %     hold on
        %     mapshow(regions);
        day
    end
end

%calcola le medie sulla 4a dimensione
if 1
    aot_dest=squeeze(nanmean(aot_dest,4));
    for i=1:365
        temp=aot_dest(:,:,i);
        temp=temp';
        aot_data(:,:,i)=temp;
    end
    lat_dest=lat_dest';
    long_dest=long_dest';
    save ..\Data\AOT_2006\aot_2006 aot_data;
    save ..\Data\AOT_2006\aot_lat lat_dest;
    save ..\Data\AOT_2006\aot_lon long_dest;
end

if 0
    for i=1:365
        for k=1:3
            if (aot_flag(i,k)==1)
                if k==1
                    figure('Position',[0,400,600,320]);
                    title(['Giorno ',num2str(i)]);
                    drawnow;
                end
                if k==2
                    figure('Position',[600,400,600,320]);
                    title(['Giorno ',num2str(i)]);
                    drawnow;
                end
                if k==3
                    figure('Position',[0,0,600,320]);
                    title(['Giorno ',num2str(i)]);
                    drawnow;
                end
                mapshow(long_dest,lat_dest,aot_dest(:,:,i,k),'DisplayType','surface');
            end
        end
        input(['Fine giorno ' num2str(i)]);
        close all;
    end
end


%filmato giorni
if 0
    for i=1:365
        figure;
        mapshow(long_dest,lat_dest,aot_data(:,:,i,1),'DisplayType','surface');
        input(['Fine giorno ' num2str(i)]);
        close all;
    end
end

