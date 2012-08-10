clc
shape=shaperead('../Maps/Scotland/scotland_only');
load('../Data/popcount_scotland2008.mat');
load('../Data/krig_lat_scotland2009');
load('../Data/krig_lon_scotland2009');
load('../Data/krig_mask_scotlandonly2009');
[lat,lon]=meshgrid(krig_lat,krig_lon);
mask(isnotnan(mask))=1;
sub=double(sub.*mask);

%plot population count
if 0
    temp=log10(sub);
    axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
        'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
    shape.X(isnan(shape.X))=0;
    shape.Y(isnan(shape.Y))=0;
    indices=find(shape.X==0);
    
    for i=1:length(indices)-1
        geoshow(shape.Y(indices(i)+1:indices(i+1)-1),shape.X(indices(i)+1:indices(i+1)-1),'DisplayType','polygon');
        hold on
    end
    
    h= geoshow(lat,lon,temp,'DisplayType','texture');
    set(h,'FaceColor','flat');
    colorbar
    %colormap('gray');
end

%plot mappa di probabilità di eccedenza e ampiezza dell'intervallo di confidenza
if 0
    alpha=0.05;
    l1=alpha/2;
    l2=1-l1;
    for t=1:1
        load('c:\Francesco\STEM3_output\multi_risk\risk_t79_level50');
        for b=1:size(thr_prob,3)
            index(t,b)=nansum(nansum(thr_prob(:,:,b).*sub));
        end
        subplot(1,2,1);
        temp=thr_prob(:,:,1);
        %temp(1,1)=0;
        %temp(1,2)=1;
        
        axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
            'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
        shape.X(isnan(shape.X))=0;
        shape.Y(isnan(shape.Y))=0;
        indices=find(shape.X==0);
        
        for i=1:length(indices)-1
            geoshow(shape.Y(indices(i)+1:indices(i+1)-1),shape.X(indices(i)+1:indices(i+1)-1),'DisplayType','polygon');
            hold on
        end
        
        h= geoshow(lat,lon,temp,'DisplayType','texture');
        set(h,'FaceColor','flat');
        colorbar
        colormap('gray');
        title(['Prob. PM_{10}>50 \mug/m^{3} - Day 1']);
        
        subplot(1,2,2);
        for j=1:size(thr_prob,2)
            for i=1:size(thr_prob,1)
                if isnotnan(thr_prob(i,j,1))
                    [f,x]=ecdf(squeeze(thr_prob(i,j,:)));
                    int1=find(f>=l1,1);
                    int2=find(f>=l2,1);
                    down=x(int1);
                    up=x(int2);
                    temp(i,j)=up-down;
                end
            end
        end

        originLat = dm2degrees([55 00]);
        originLon = dm2degrees([-4 0]);
        axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
            'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
        shape.X(isnan(shape.X))=0;
        shape.Y(isnan(shape.Y))=0;
        indices=find(shape.X==0);
        
        for i=1:length(indices)-1
            geoshow(shape.Y(indices(i)+1:indices(i+1)-1),shape.X(indices(i)+1:indices(i+1)-1),'DisplayType','polygon');
            hold on
        end
        
        h= geoshow(lat,lon,temp,'DisplayType','texture');
        set(h,'FaceColor','flat');
        colorbar
        colormap('gray');
        %title(['Population count distribution']);
    end
end

%plot dell'indice di rischio giornaliero e intervallo di confidenza
if 0
    alpha=0.05;
    l1=alpha/2;
    l2=1-l1;
    files=dir('c:\Francesco\STEM3_output\risk\PM10\*.mat');
    for i=1:length(files)
        [i length(files)]
        name=files(i).name;
        load(['c:\Francesco\STEM3_output\risk\PM10\',name]);
        index=[];
        for b=1:size(thr_prob,3)
            index(b)=nansum(nansum(thr_prob(:,:,b).*sub));
        end
        indices=strfind(name,'_');
        day=str2num(name(indices(1)+2:indices(2)-1));
        [f,x]=ecdf(index);
        media(day)=mean(index);
        int1=find(f>=l1,1);
        int2=find(f>=l2,1);
        down(day)=x(int1);
        up(day)=x(int2);
    end
    plot(log((up+down)/2));
    hold on
    plot(log(down),':');
    plot(log(up),':');
    %errorbar(log10((up+down)/2),(log10(up)-log10(down))/2,'x-');
    errorbar((up+down)/2,(up-down)/2,'x');
end

% figure
% hold on
% for i=1:length(down)
%     x1=i;
%     x2=i;
%     y1=down(i);
%     y2=up(i);
%     plot([x1,x2],[y1,y2],'-');
%     plot(i,(y1+y2)/2,'.');
% end

%plot probabilità di superamento giorni di eccedenza
if 0
    files=dir('c:\Francesco\STEM3_output\risk\*.mat');
    for i=1:length(files)
        [i length(files)]
        name=files(i).name;
        load(['c:\Francesco\STEM3_output\risk\',name]);
        indices=strfind(name,'_');
        day=str2num(name(indices(1)+2:indices(2)-1));
        prob(:,:,day)=mean(thr_prob,3);
    end
end

if 0
    load prob
    nsim=1000;
    nday=zeros(size(prob,1),size(prob,2),nsim);
    for b=1:nsim %level 50 for pm10
        b
        temp=zeros(size(prob,1),size(prob,2),size(prob,3));
        for t=1:size(prob,3)
            temp(:,:,t)=binornd(1,prob(:,:,t));
        end
        nday(:,:,b)=sum(temp,3);
    end
    save nday nday -v7.3
end

if 0
    prob=zeros(size(nday,1),size(nday,2));
    for j=1:size(nday,2)
        j
        for i=1:size(nday,1)
            data=squeeze(nday(i,j,:));
            L=data>7;
            prob(i,j)=sum(L)/length(L);
        end
    end
end
figure
if 1
    subplot(1,2,1);
    temp=diff;
    %temp(1,1)=0;
    %temp(1,2)=1;
    
    axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
        'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
    shape.X(isnan(shape.X))=0;
    shape.Y(isnan(shape.Y))=0;
    indices=find(shape.X==0);
    
    for i=1:length(indices)-1
        geoshow(shape.Y(indices(i)+1:indices(i+1)-1),shape.X(indices(i)+1:indices(i+1)-1),'DisplayType','polygon');
        hold on
    end
    
    h= geoshow(lat,lon,temp,'DisplayType','texture');
    set(h,'FaceColor','flat');
    colorbar
    %colormap('gray');
    %title(['Prob. PM_{10}>50 \mug/m^{3} - Day 1']);
    
    subplot(1,2,2);
    temp=Var_PM25_multi;
    
    originLat = dm2degrees([55 00]);
    originLon = dm2degrees([-4 0]);
    axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
        'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
    shape.X(isnan(shape.X))=0;
    shape.Y(isnan(shape.Y))=0;
    indices=find(shape.X==0);
    
    for i=1:length(indices)-1
        geoshow(shape.Y(indices(i)+1:indices(i+1)-1),shape.X(indices(i)+1:indices(i+1)-1),'DisplayType','polygon');
        hold on
    end
    
    h= geoshow(lat,lon,temp,'DisplayType','texture');
    set(h,'FaceColor','flat');
    colorbar
    %colormap('gray');
    %title(['Population count distribution']);
end
