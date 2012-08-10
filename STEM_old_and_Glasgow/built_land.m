clc
clear all
load land_use_all_layers

if 1
    lat_start=43;
    lat_end=47;
    lon_start=6;
    lon_end=14;
    step=0.05;
    lat=[lat_start:step:lat_end];
    lon=[lon_start:step:lon_end];
    land_use=zeros(length(lon),length(lat),365);
    %leggo i nomi dei file nella cartella dei dati
    if 0
        cartella='C:\Documents and Settings\Francesco\Desktop\Land\*.hdf';
        files=dir(cartella);
        for i=1:length(files)
            fname=[cartella(1:end-5) files(i).name];
            data = hdfread(fname,'MOD12C1','Fields','Land_Cover_Type_1_Percent');
            lat_index=round((90-lat_end)/0.05):1:round((90-lat_start)/0.05);
            lon_index=round((180+lon_start)/0.05):1:round((180+lon_end)/0.05);
            data=double(data(lon_index,lat_index,14)/255*100)
            land_use=repmat(data,[1 1 365]);
        end
    end
    land_use=double(repmat(land_use_all_layers(:,:,14),[1 1 365]));
    save land_use land_use
    save lat_land lat
    save lon_land lon
end