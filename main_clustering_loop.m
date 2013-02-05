%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc
%clear all
%close all

%RandStream.setDefaultStream(RandStream('mt19937ar','seed',2222));

flag_lakedata=0;

flag_parallel=0;
flag_covariates=0;
flag_remove_zero=0;
flag_simulatedata=0;
flag_space=0;
pathparallel='/opt/matNfs/';

max_log=-999999;
n_clusters=8;

for zz=1:50
    if flag_simulatedata
        %Data simulation
        phase_step=10; %degree
        t=1:260;
        if 1
            clusters=[80 40 20 10 5];
            noise_std=[0.5,0.4,0.3,0.2,0.1];
            phase=0;
            Y=[];
            counter=1;
            for i=1:length(clusters)
                for j=1:clusters(i)
                    Y(counter,:)=sin(t/8.276+phase/360*2*pi)+normrnd(0,noise_std(i),1,length(t));
                    counter=counter+1;
                end
                phase=phase+phase_step;
            end
        else
            counter=1;
            for i=1:20
                Y(counter,:)=sin(t/12)+normrnd(0,noise_std,1,length(t))+t*0.005;
                counter=counter+1;
            end
            for i=1:20
                Y(counter,:)=sin(t/12+phase_step/360*2*pi)+normrnd(0,noise_std,1,length(t))+t*0.005;
                counter=counter+1;
            end
            for i=1:20
                Y(counter,:)=sin(t/12)+normrnd(0,noise_std,1,length(t))-t*0.005;
                counter=counter+1;
            end
            for i=1:20
                Y(counter,:)=sin(t/12+phase_step/360*2*pi)+normrnd(0,noise_std,1,length(t))-t*0.005;
                counter=counter+1;
            end
        end
        sd_g.Y{1}=Y;
        sd_g.Y_name={'random_data'};
        N=size(Y,1);
        %it is not important where the sites are so their coordinates are randomly generated
        st_grid=stem_grid(rand(N,2),'deg','sparse','point');
        shape=[];
    else
        if flag_lakedata
            load C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\arclake\lake_average\lakes_data
            load C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\arclake\lake_average\lakes_covariates.mat
            idx=[];
            id_vec=cell2mat(lakes_covariate(:,1));
            for i=1:length(data.id_day)
                idx=[idx;find(id_vec==data.id_day(i))];
            end
            X=cell2mat(lakes_covariate(idx,[5,6,7,14]));
            X_label=lakes_covariate_label([5,6,7,14]);
            
            shape = shaperead('C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\arclake\lake_average\CNTR_BN_03M_2010');
            temp=data.day;
            if flag_remove_zero
                temp(temp<=273.16)=NaN;
            end
            load C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\arclake\lake_average\residual_timeseries
            sd_g.Y{1}=temp;
            %sd_g.Y{1}=diff_all;
            sd_g.Y_name={'temperature'};
            %st_grid=stem_grid(data.coordinates_day(1:end-1,:),'deg','sparse','point');
            st_grid=stem_grid(data.coordinates_day,'deg','sparse','point');
        else
            load C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\TOS\TOS_dataset
            shape = shaperead('C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\TOS\Maps\scotland_only');
            shape_water = shaperead('C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\TOS\Maps\waterways');
            shape_natural = shaperead('C:\Francesco\Universita\Ricerca\Visiting\2012_Glasgow\ricerca\TOS\Maps\natural');
            Y=data.Y_rivers_monthly_small;
            coordinates=data.river_coordinates;
            sd_g.Y{1}=Y;
            sd_g.Y_name={'TOC'};
            st_grid=stem_grid(coordinates,'deg','sparse','point');
        end
    end
    
    %Model definition

    N=size(sd_g.Y{1},1);
    T=size(sd_g.Y{1},2);
    w=unifrnd(0,1,N,n_clusters);
    s=sum(w,2);
    sd_g.X_time_name{1}=[];
    for i=1:n_clusters
        w(:,i)=w(:,i)./s;
        sd_g.X_time_name{1}{i}={['weights_',num2str(i)]};
    end
    sd_g.X_time{1}=w;
    
    if flag_covariates
        sd_g.X_beta{1}=X;
        sd_g.X_beta_name{1}=X_label;
    else
        sd_g.X_beta=[];
        sd_g.X_beta_name=[];
    end
    
    if flag_space
        sd_g.X_g{1}=ones(N,1,1,1);
        sd_g.X_g_name{1}={'constant'};
    else
        sd_g.X_g=[];
        sd_g.X_g_name=[];
    end
    
    st_varset_g=stem_varset(sd_g.Y,sd_g.Y_name,[],[],sd_g.X_beta,sd_g.X_beta_name,sd_g.X_time,sd_g.X_time_name,sd_g.X_g,sd_g.X_g_name);
    st_gridlist_g=stem_gridlist();
    st_gridlist_g.add(st_grid);
    
    %datestamp
    if not(flag_simulatedata)
        if flag_lakedata
            st_datestamp=stem_datestamp(data.time(1),data.time(end),T);
        else
            st_datestamp=stem_datestamp(1,T,T);
        end
        st_data=stem_data(st_varset_g,st_gridlist_g,[],[],st_datestamp,shape);
        if flag_lakedata
            st_data.time_crop(T-365*5:T);
            %st_data.time_crop(1:365*5);
            st_data.time_average(7);
        end
        st_data.standardize_sbs;
        %st_data.detrend;
    else
        st_datestamp=stem_datestamp(1,T,T);
        st_data=stem_data(st_varset_g,st_gridlist_g,[],[],st_datestamp,shape);
    end
    
    time_diagonal=1;
    clustering=1;
    theta_clustering=10;
    st_par=stem_par(st_data,'exponential',0,time_diagonal,clustering,theta_clustering);
    st_model=stem_model(st_data,st_par);
    st_model.set_system_size();
    
    %st_par initialization
    if flag_space
        st_par.alpha_g=0.5;
        st_par.theta_g=100;
        st_par.v_g=1;
    end
    st_par.sigma_eta=diag(repmat(0.5,n_clusters,1)+unifrnd(-0.4,0.4,n_clusters,1));
    st_par.G=diag(repmat(0.8,n_clusters,1)+unifrnd(-0.2,0.2,n_clusters,1));
    st_par.sigma_eps=diag(0.2)+unifrnd(-0.1,0.3,1,1);
    st_model.set_initial_values(st_par);
    
    %model estimation
    st_EM_options=stem_EM_options(0.001,30,'single',[],0);
    st_model.EM_estimate(st_EM_options);
    st_model.set_logL;
    ll=st_model.stem_EM_result.logL
    cc=round(sum(st_model.stem_data.X_time))
    if ll>max_log
        save st_model_max2 st_model
        max_log=ll;
    end
    [zz max_log]
end

%result evaluation
if 0
    [val,idx]=max(st_model.stem_data.X_time,[],2);
    if not(flag_simulatedata)
        figure
        if flag_lakedata
            mapshow(st_model.stem_data.shape,'Color', 'black');
        else
            axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
                'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
            st_model.stem_data.shape.X(isnan(st_model.stem_data.shape.X))=0;
            st_model.stem_data.shape.Y(isnan(st_model.stem_data.shape.Y))=0;
            indices=find(st_model.stem_data.shape.X==0);
            for j=1:length(indices)-1
                geoshow(st_model.stem_data.shape.Y(indices(j)+1:indices(j+1)-1),st_model.stem_data.shape.X(indices(j)+1:indices(j+1)-1),'DisplayType','polygon','FaceColor','none');
            end
            
            %     for j=1:length(shape_water)
            %         geoshow(shape_natural(j).Y,shape_natural(j).X,'DisplayType','polygon','FaceColor','none','EdgeColor','blue');
            %     end
        end
        hold on
        
        for i=1:length(idx)
            switch idx(i)
                case 1
                    color='r';
                    mark='o';
                case 2
                    color='g';
                    mark='o';
                case 3
                    color='b';
                    mark='o';
                case 4
                    color='c';
                    mark='o';
                case 5
                    color='m';
                    mark='o';
                case 6
                    color='y';
                    mark='o';
                case 7
                    color='y';
                    mark='^';
                case 8
                    color='r';
                    mark='^';
                case 9
                    color='b';
                    mark='^';
                case 10
                    color='g';
                    mark='^';
                case 11
                    color='m';
                    mark='^';
                case 12
                    color='c';
                    mark='^';
            end
            if flag_lakedata
                msize=round(val(i)*15);
                plot(st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(i,2),st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(i,1),[color,mark],'MarkerFaceColor',color,'MarkerEdgeColor','k','MarkerSize',10);
            else
                msize=round(val(i)*10);
                geoshow(st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(i,1),st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(i,2),'DisplayType','Point','Marker','o','MarkerFaceColor',color,'MarkerEdgeColor','k','MarkerSize',msize);
            end
        end
        % [lon_click,lat_click] = ginput(1);
        % d=(st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(:,1)-lat_click).^2+(st_model.stem_data.stem_gridlist_g.grid{1}.coordinate(:,2)-lon_click).^2;
        % [val2,idx2]=min(d);
        % figure
        % plot(st_model.stem_EM_result.stem_kalmansmoother_result.zk_s(idx(idx2),2:end)','LineWidth',2);
        % hold on
        % plot(st_model.stem_data.Y(idx2,:),'r','LineWidth',2);
    end
    
    figure
    H=ceil(sqrt(n_clusters));
    K=ceil(n_clusters/H);
    counter=1;
    for h=1:H
        for k=1:K
            if counter<=n_clusters
                subplot(H,K,counter);
                hold on
                counter2=0;
                mult_value=0;
                for i=1:length(idx)
                    if idx(i)==counter
                        plot(st_model.stem_data.Y(i,:),'LineWidth',1,'Color','r');
                        counter2=counter2+1;
                        mult_value=mult_value+abs(st_model.stem_data.X_time(i,idx(i)));
                    end
                end
                mult_value=mult_value/counter2;
                mult_value=1;
                plot(mult_value*st_model.stem_EM_result.stem_kalmansmoother_result.zk_s(counter,2:end)','LineWidth',3,'Color','b');
                temp1=mult_value*(st_model.stem_EM_result.stem_kalmansmoother_result.zk_s(counter,2:end)'+3*sqrt(squeeze(st_model.stem_EM_result.stem_kalmansmoother_result.Pk_s(counter,counter,2:end))));
                temp2=mult_value*(st_model.stem_EM_result.stem_kalmansmoother_result.zk_s(counter,2:end)'-3*sqrt(squeeze(st_model.stem_EM_result.stem_kalmansmoother_result.Pk_s(counter,counter,2:end))));
                plot(temp1,'LineWidth',2,'Color','b','LineStyle','--');
                plot(temp2,'LineWidth',2,'Color','b','LineStyle','--');
                title(['Cluster ',num2str(counter),' - ',num2str(counter2),' time series']);
                ylim([-5,5]);
                grid on
            end
            counter=counter+1;
        end
    end
    set(gcf, 'Position', get(0,'Screensize'));
    
    %     figure
    %     counter=1;
    %     diff_all=[];
    %     for h=1:H
    %         for k=1:K
    %             if counter<=n_clusters
    %                 cluster{counter}=[];
    %                 subplot(H,K,counter);
    %                 hold on
    %                 counter2=0;
    %                 for i=1:length(idx)
    %                     if idx(i)==counter
    %                         diff=st_model.stem_data.Y(i,:)-st_model.stem_EM_result.stem_kalmansmoother_result.zk_s(counter,2:end);
    %                         cluster{counter}=[cluster{counter};diff(:)];
    %                         diff_all=[diff_all;diff];
    %                         plot(diff,'LineWidth',1,'Color','r');
    %                         counter2=counter2+1;
    %                     end
    %                 end
    %                 title([num2str(counter2),' time series in this cluster']);
    %                 ylim([-5,5]);
    %                 grid on
    %             end
    %             counter=counter+1;
    %         end
    %     end
    %     set(gcf, 'Position', get(0,'Screensize'));
    %
    %     figure
    %     counter=1;
    %     for h=1:H
    %         for k=1:K
    %             if counter<=n_clusters
    %                 subplot(H,K,counter);
    %                 hold on
    %                 if not(isempty(cluster{counter}))
    %                     histfit(cluster{counter},30);
    %                 end
    %             end
    %             counter=counter+1;
    %         end
    %     end
    %     set(gcf, 'Position', get(0,'Screensize'));
end





