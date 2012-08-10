classdef stem_bootkrig_result < handle
    
    properties
        stem_grid=[];
        variable_name=[];
        Y_hat=[];
        Var_Y_hat=[];
        shape=[];
    end
    
    methods
        function obj = stem_bootkrig_result(variable_name,stem_grid,shape)
            if nargin<3
                error('Not enough input arguments');
            end
            obj.variable_name=variable_name;
            obj.stem_grid=stem_grid;
            obj.shape=shape;
        end
        
        function plot(obj,time_step)
            if nargin<1
                error('Not enough input arguments');
            end
            if (time_step<1)||(time_step>size(obj.Y_hat,3))
                error('time_step out of bounds');
            end

            figure
            for t=time_step
                data_mean=nanmean(obj.Y_hat(:,:,t,:),4);
                data_var=nanvar(obj.Y_hat(:,:,t,:),[],4);
                lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                long=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                subplot(1,2,1);
                mapshow(obj.shape);
                hold on
                %colormap('gray');
                mapshow(long,lat,data_mean,'DisplayType','texture');
                set(gca,'Xlim',[obj.stem_grid.box(3),obj.stem_grid.box(4)]);
                set(gca,'Ylim',[obj.stem_grid.box(1),obj.stem_grid.box(2)]);
                colorbar
                title(['Bootstrap average kriged ',obj.variable_name,' - day ',num2str(t)]);
                subplot(1,2,2);
                mapshow(obj.stem_model.stem_data.shape);
                hold on
                %colormap('gray');
                mapshow(long,lat,data_var,'DisplayType','texture');
                set(gca,'Xlim',[obj.stem_grid.box(3),obj.stem_grid.box(4)]);
                set(gca,'Ylim',[obj.stem_grid.box(1),obj.stem_grid.box(2)]);
                colorbar
                title(['Bootstrap average kriging variance ',obj.variable_name,' - day ',num2str(t)]);
                pause(0.5);
            end
        end
        
        function map = bootstrap_prob(obj,time_step,type,value,graph)
            %Plot for a particular time step the map probability that the
            %variable_name variable is either lower or higher a value
            % P(variable_name(s,t)><value)
            
            % variable_name:    variable to be plotted
            % time_step:        time step
            % type:             either 'lower' or 'higher'
            % value:            any number
            if nargin<4
                error('Not enough input parameters');
            end
            if nargin<5
                graph=0;
            end
            if (time_step<1)||(time_step>size(obj.Y_hat,3))
                error('time_step out of bounds');
            end
            if sum(strcmp(type,{'lower','higher'}))==0
                error('The type can only be lower or higher');
            end
            if isscalar(value)==0
                error('value must be a scalar');
            end
            
            map=zeros(size(obj.Y_hat,1),size(obj.Y_hat,2));
            for j=1:size(obj.Y_hat,2)
                for i=1:size(obj.Y_hat,1)
                    data=squeeze(obj.Y_hat(i,j,time_step,:));
                    if sum(isnotnan(data))>0
                        if strcmp(type,'lower')
                            if min(data)>value
                                map(i,j)=0;
                            else
                                if max(data)<value
                                    map(i,j)=1;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-value));
                                    map(i,j)=x(index);
                                end
                            end
                        else
                            if min(data)>value
                                map(i,j)=1;
                            else
                                if max(data)<value
                                    map(i,j)=0;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-value));
                                    map(i,j)=1-x(index);
                                end
                            end
                        end
                    else
                        map(i,j)=NaN;
                    end
                end
            end
            if graph
                lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                long=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                mapshow(obj.shape);
                hold on
                %colormap('gray');
                mapshow(long,lat,map,'DisplayType','texture');
                set(gca,'Xlim',[obj.stem_grid.box(3),obj.stem_grid.box(4)]);
                set(gca,'Ylim',[obj.stem_grid.box(1),obj.stem_grid.box(2)]);
                xlabel('Longitude');
                ylabel('Latitude');
                if strcmp(type,'higher')
                    simble=' > ';
                else
                    simble=' < ';
                end
                title(['Prob. ',obj.variable_name,simble,num2str(value),' day ',num2str(t)]);
                colorbar
            end
        end
        
        function [map, nstepdist] = bootstrap_prob_timestep(obj,type,value,thr_type,thr,graph)
            %Plot the mean number of time steps for which the
            %variable_name variable is above/below a value
            %Plot the map probability that the number of time steps for
            %which variable_name variable is above/below a value is
            %above/below a threshold thr.
            %variable_name:     the variable to be plotted
            %type:              either 'lower' or 'higher'
            %value:
            %thr_type:          either 'lower' or 'higher'
            %thr:               number of time steps
            if nargin<5
                error('Not enough input parameters');
            end
            if nargin<6
                graph=0;
            end
            if sum(strcmp(type,{'lower','higher'}))==0
                error('The type can only be lower or higher');
            end
            if isscalar(value)==0
                error('value must be a scalar');
            end
            if sum(strcmp(thr_type,{'lower','higher'}))==0
                error('The thr_type can only be lower or higher');
            end
            if (thr<1)||(thr>size(obj.Y_hat,3))
                error('thr cannot be lower than 1 and higher than T');
            end
            if thr~=round(thr)
                error('thr must be an integer value');
            end
            map=zeros(size(obj.Y_hat,1),size(obj.Y_hat,2));
            nstepdist=zeros(size(obj.Y_hat,1),size(obj.Y_hat,2),size(obj.Y_hat,4));
            
            for b=1:size(obj.Y_hat,4)
                data=zeros(size(obj.Y_hat,1),size(obj.Y_hat,2));
                for t=1:size(obj.Y_hat,3)
                    if strcmp(type,'lower')
                        data=data+(obj.Y_hat(:,:,t,b)<value);
                    else
                        data=data+(obj.Y_hat(:,:,t,b)>value);
                    end
                end
                %nstepdist is composed of b maps. Along the 3rd dimension
                %is the distribution of the number of days the pixel
                %was above (or below) the threshold.
                nstepdist(:,:,b)=data;
            end
            
            for j=1:size(obj.Y_hat,2)
                for i=1:size(obj.Y_hat,1)
                    data=squeeze(nstepdist(i,j,:));
                    if sum(isnotnan(data))>0
                        if strcmp(thr_type,'lower')
                            if min(data)>thr
                                map(i,j)=0;
                            else
                                if max(data)<thr
                                    map(i,j)=1;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-thr));
                                    map(i,j)=x(index);
                                end
                            end
                        else
                            if min(data)>thr
                                map(i,j)=1;
                            else
                                if max(data)<thr
                                    map(i,j)=0;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-thr));
                                    map(i,j)=1-x(index);
                                end
                            end
                        end
                    else
                        map(i,j)=NaN;
                    end
                end
            end
            
            if graph
                lat=reshape(obj.stem_grid.coordinate(:,1),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                long=reshape(obj.stem_grid.coordinate(:,2),obj.stem_grid.grid_size(1),obj.stem_grid.grid_size(2));
                subplot(1,2,1);
                mapshow(obj.shape);
                hold on
                mapshow(long,lat,nanmean(nstepdist,3),'DisplayType','texture');
                set(gca,'Xlim',[obj.stem_grid.box(3),obj.stem_grid.box(4)]);
                set(gca,'Ylim',[obj.stem_grid.box(1),obj.stem_grid.box(2)]);
                xlabel('Longitude');
                ylabel('Latitude');
                if strcmp(type,'higher')
                    simble=' > ';
                else
                    simble=' < ';
                end
                title(['Mean # of time steps for which ',obj.variable_name,simble,num2str(value)]);
                colorbar
                subplot(1,2,2);
                mapshow(obj.shape);
                hold on
                %colormap('gray');
                mapshow(long,lat,map,'DisplayType','texture');
                set(gca,'Xlim',[obj.stem_grid.box(3),obj.stem_grid.box(4)]);
                set(gca,'Ylim',[obj.stem_grid.box(1),obj.stem_grid.box(2)]);
                xlabel('Longitude');
                ylabel('Latitude');
                if strcmp(type,'higher')
                    simble=' > ';
                else
                    simble=' < ';
                end
                if strcmp(thr_type,'higher')
                    thr_simble=' > ';
                else
                    thr_simble=' < ';
                end
                title(['Prob. # time steps for which ',obj.variable_name,simble,num2str(value),thr_simble,num2str(thr)]);
                colorbar
            end
        end
    end
    
    methods (Static)
        function [out,residuals] = ma(in,window)
            if nargin<2
                window=10;
            end
            if not(isvector(in))
                error('The data input parameter must be a vector');
            end
            out=zeros(1,length(in));
            for i=1:length(in)
                idx1=i-window;
                idx2=i+window;
                if idx1<1
                    idx1=1;
                end
                if idx2>length(in)
                    idx2=length(in);
                end
                out(i)=nanmean(in(idx1:idx2));
            end
            residuals=in-out;
        end
        
        function map = bootstrap_prob_lff(directory,variable_name,time_step,type,value,graph,shape)
            %Plot for a particular time step the map probability that the
            %variable_name variable is either lower or higher a value
            % P(variable_name(s,t)><value)
            
            % variable_name:    variable to be plotted
            % time_step:        time step
            % type:             either 'lower' or 'higher'
            % value:            any number
            if nargin<5
                error('Not enough input parameters');
            end
            if nargin<6
                graph=0;
            end
            if (graph==1)&&(nargin==6)
                error('A shape file must be provided for the plot');
            end
            if sum(strcmp(type,{'lower','higher'}))==0
                error('The type can only be lower or higher');
            end
            if isscalar(value)==0
                error('value must be a scalar');
            end
            if time_step<1
                error('time_step must be > 0');
            end
            if not(strcmp(directory(end),'/'))
                directory=[directory,'/'];
            end
            %reading the file containing grid information
            files=dir([directory,'t1/*.mat']);
            if isempty(files)
                error('The directory for t=1 cannot be found. Grid information cannot be extracted');
            end
            load([directory,'t1/',files(1).name]);
            Grid=boot_krig_fixtime.Grid;
            Grid_size=boot_krig_fixtime.Grid_size;
            if not(strcmp(boot_krig_fixtime.variable_name,variable_name))
                error('Files refers to a different varible or the variable name provided is incorrect');
            end
            files=dir([directory,'t',num2str(time_step),'/*.mat']);
            if isempty(files)
                error(['The directory for time step ',num2str(time_step),' is empty or the directory name is incorrect']);
            end
            for b=1:length(files)
                disp(['Loading ',num2str(b),' of ',num2str(length(files))]);
                load([directory,'t',num2str(time_step),'/',files(b).name]);
                boot_data(:,:,b)=boot_krig_fixtime.data;
            end
            
            map=zeros(size(boot_data,1),size(boot_data,2));
            
            for j=1:size(boot_data,2)
                for i=1:size(boot_data,1)
                    data=squeeze(boot_data(i,j,:));
                    if sum(isnotnan(data))>0
                        if strcmp(type,'lower')
                            if min(data)>value
                                map(i,j)=0;
                            else
                                if max(data)<value
                                    map(i,j)=1;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-value));
                                    map(i,j)=x(index);
                                end
                            end
                        else
                            if min(data)>value
                                map(i,j)=1;
                            else
                                if max(data)<value
                                    map(i,j)=0;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-value));
                                    map(i,j)=1-x(index);
                                end
                            end
                        end
                    else
                        map(i,j)=NaN;
                    end
                end
            end
            if graph
                lat=reshape(Grid(:,1),Grid_size(2),Grid_size(1));
                long=reshape(Grid(:,2),Grid_size(2),Grid_size(1));
                mapshow(shape);
                hold on
                %colormap('gray');
                mapshow(long,lat,map,'DisplayType','texture');
                set(gca,'Xlim',[min(long(:)),max(long(:))]);
                set(gca,'Ylim',[min(lat(:)),max(lat(:))]);
                xlabel('Longitude');
                ylabel('Latitude');
                if strcmp(type,'higher')
                    simble=' > ';
                else
                    simble=' < ';
                end
                title(['Prob. ',variable_name,simble,num2str(value),' day ',num2str(time_step)]);
                colorbar
            end
        end
        
        function [map, nstepdist] = bootstrap_prob_timestep_lff(directory,variable_name,type,value,thr_type,thr,graph,shape)
            %Plot the mean number of time steps for which the
            %variable_name variable is above/below a value
            %Plot the map probability that the number of time steps for
            %which variable_name variable is above/below a value is
            %above/below a threshold thr.
            %variable_name:     the variable to be plotted
            %type:              either 'lower' or 'higher'
            %value:
            %thr_type:          either 'lower' or 'higher'
            %thr:               number of time steps
            if nargin<6
                error('Not enough input parameters');
            end
            if nargin<7
                graph=0;
            end
            if (graph==1)&&(nargin==7)
                error('A shape file must be provided for the plot');
            end
            if sum(strcmp(type,{'lower','higher'}))==0
                error('The type can only be lower or higher');
            end
            if isscalar(value)==0
                error('value must be a scalar');
            end
            if sum(strcmp(thr_type,{'lower','higher'}))==0
                error('The thr_type can only be lower or higher');
            end
            if (thr<1)
                error('thr cannot be lower than 1');
            end
            if thr~=round(thr)
                error('thr must be an integer value');
            end
            if not(strcmp(directory(end),'/'))
                directory=[directory,'/'];
            end
            %count the number of timesteps
            files=dir(directory);
            T=0;
            for i=1:length(files)
                if (files(i).isdir)&&(strcmp(files(i).name(1),'t'))
                    T=T+1;
                end
            end
            if T==0
                error('The directory does not contain any valid sub-directory');
            end
            
            %reading the file containing grid information and count the
            %bootstrap recplications
            files=dir([directory,'t1/*.mat']);
            if isempty(files)
                error('The directory for t=1 cannot be found. Grid information cannot be extracted');
            end
            B=length(files);
            load([directory,'t1/',files(1).name]);
            Grid=boot_krig_fixtime.Grid;
            Grid_size=boot_krig_fixtime.Grid_size;
            if not(strcmp(boot_krig_fixtime.variable_name,variable_name))
                error('Files refers to a different varible or the variable name provided is incorrect');
            end
            
            map=zeros(Grid_size(2),Grid_size(1));
            nstepdist=zeros(Grid_size(2),Grid_size(1),B);
            
            for b=1:B
                data=zeros(Grid_size(2),Grid_size(1));
                for t=1:T
                    disp(['Loading file t ',num2str(t),'/',num2str(T),' - b ',num2str(b),'/',num2str(B)]);
                    load([directory,'t',num2str(t),'/',files(b).name]);
                    if not(boot_krig_fixtime.temporal_index==t)
                        error('Temporal index differs from what expected');
                    end
                    if not(boot_krig_fixtime.bootstrap_replica_index==b)
                        error('Bootstrap replica index differs from what expected');
                    end
                    if strcmp(type,'lower')
                        data=data+(boot_krig_fixtime.data<value);
                    else
                        data=data+(boot_krig_fixtime.data>value);
                    end
                end
                %nstepdist is composed of b maps. Along the 3rd dimension
                %is the distribution of the number of days the pixel
                %was above (or below) the threshold.
                nstepdist(:,:,b)=data;
            end
            
            disp('Evaluating distributions...');
            for j=1:Grid_size(1)
                for i=1:Grid_size(2)
                    data=squeeze(nstepdist(i,j,:));
                    if sum(isnotnan(data))>0
                        if strcmp(thr_type,'lower')
                            if min(data)>thr
                                map(i,j)=0;
                            else
                                if max(data)<thr
                                    map(i,j)=1;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-thr));
                                    map(i,j)=x(index);
                                end
                            end
                        else
                            if min(data)>thr
                                map(i,j)=1;
                            else
                                if max(data)<thr
                                    map(i,j)=0;
                                else
                                    [x,f]=ecdf(data);
                                    [val,index]=min(abs(f-thr));
                                    map(i,j)=1-x(index);
                                end
                            end
                        end
                    else
                        map(i,j)=NaN;
                    end
                end
            end
            disp('Evaluation ended');
            
            if graph
                %                 lat=reshape(obj.Grid(:,1),obj.Grid_size(1),obj.Grid_size(2));
                %                 long=reshape(obj.Grid(:,2),obj.Grid_size(1),obj.Grid_size(2));
                %                 subplot(1,2,1);
                mapshow(shape);
                %                 hold on
                %                 mapshow(long,lat,nanmean(nstepdist,3),'DisplayType','texture');
                %                 set(gca,'Xlim',[obj.box(3),obj.box(4)]);
                %                 set(gca,'Ylim',[obj.box(1),obj.box(2)]);
                %                 xlabel('Longitude');
                %                 ylabel('Latitude');
                %                 if strcmp(type,'higher')
                %                     simble=' > ';
                %                 else
                %                     simble=' < ';
                %                 end
                %                 title(['Mean # of time steps for which ',variable_name,simble,num2str(value)]);
                %                 colorbar
                %                 subplot(1,2,2);
                %                 mapshow(obj.stem_model.stem_data.shape);
                %                 hold on
                %                 %colormap('gray');
                %                 mapshow(long,lat,map,'DisplayType','texture');
                %                 set(gca,'Xlim',[obj.box(3),obj.box(4)]);
                %                 set(gca,'Ylim',[obj.box(1),obj.box(2)]);
                %                 xlabel('Longitude');
                %                 ylabel('Latitude');
                %                 if strcmp(type,'higher')
                %                     simble=' > ';
                %                 else
                %                     simble=' < ';
                %                 end
                %                 if strcmp(thr_type,'higher')
                %                     thr_simble=' > ';
                %                 else
                %                     thr_simble=' < ';
                %                 end
                %                 title(['Prob. # time steps for which ',variable_name,simble,num2str(value),thr_simble,num2str(thr)]);
                %                 colorbar
            end
        end

    end
end