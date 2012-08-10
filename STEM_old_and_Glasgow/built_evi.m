clc
clear all
%carica la mappa di destinazione

if 1
    lat_start=44;
    lat_end=47;
    lon_start=6;
    lon_end=14;
    step=0.05;
    lat=[lat_start:step:lat_end];
    lon=[lon_start:step:lon_end];
    evi=zeros(length(lat),length(lon),365);
    %leggo i nomi dei file nella cartella dei dati
    cartella='C:\Users\Francesco\Desktop\NDVI\OK\*.hdf';
    files=dir(cartella);
    evi_dest=zeros(length(lat),length(lon),365)*NaN;
    day_index=1;
    counter=0;
    for i=1:length(files)
        day_start=str2num(files(i).name(14:16));
        day_end=day_start+16-1;
        fname=[cartella(1:end-5) files(i).name];
        data = hdfread(fname,'MODIS_Grid_16Day_VI_CMG','Fields','CMG 0.05
Deg 16 days EVI');
        lat_index=round((90-lat_end)/0.05):1:round((90-lat_start)/0.05);
        lon_index=round((180+lon_start)/0.05):1:round((180+lon_end)/0.05);
        data=double(data(lat_index,lon_index));
        evi(:,:,day_start:day_end)=repmat(data,[1 1 16]);
        i
    end
end
