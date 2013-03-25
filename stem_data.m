%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_data < handle
    
    %CONSTANTS
    %N  = n1_g+...+nq_g+n1_r+...+nq_r - total number of observation sites
    %Ng = n1_g+...+nq_g - total number of ground level sites
    %Nr = n1_r+...+nq_r - total number of remote sensing sites
    %M = 2 if both ground level and remote sense data are considered. M = 1 if only ground level data are considered.
    %T - number of temporal steps
    %TT = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    
    properties
        stem_varset_g=[];       %[stem_varset object] (1x1) stem_varset object for the ground level variables
        stem_varset_r=[];       %[stem_varset object] (1x1) stem_varset object for the remote sensing variables
        stem_gridlist_g=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the ground level variables        
        stem_gridlist_r=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the remote sensing variables
        stem_datestamp=[]       %[stem_datestamp object] (1x1) stem_datestamp object with information on time steps
        stem_crossval=[];       %[stem_crossval object] (1x1) stem_crossval object with information on crossvalidation
        
        shape=[];               %[struct] geographic data structure loaded from a shapefile with the boundary of the geographic region
        simulated=0;            %[boolean] (1x1) 1 if the data have been simulated, 0 otherwise.
        remote_correlated=0;    %[boolean] (1x1) 1 if the remote variables are cross-correlated
        X_time=[];              %[double]     (NxpxTT) the full X_time matrix
    end
    
    properties (SetAccess = private) 
        Y=[];                   %[double]     (NxT) the full observation matrix
        X_rg=[];                %[double]     (Nx1XTT) the full X_rg matrix
        X_beta=[];              %[double]     (NxBxTT) the full X_beta matrix
        X_g=[];                 %[double]     (Nx1xTTxK) the full X_g matrix
        DistMat_g=[];           %[double]     (Ng x Ng) distance matrix of the ground level sites
        DistMat_r=[];           %[double]     (Nr x Nr) distance matrix of the remote sensing sites
        M=[];                   %[integer >1] (Ngx1) vector of indices of the remote sensing pixels related to the ground level sites
        %flags
        can_reset=0;            %[boolean]    (1x1) 1 if data are saved to disk and can be reloaded
        X_rg_tv=0;              %[boolean]    (1x1) 1 if X_rg is time variant 0 otherwise
        X_beta_tv=0;            %[boolean]    (1x1) 1 if X_beta is time variant 0 otherwise
        X_time_tv=0;            %[boolean]    (1x1) 1 if X_time is time variant 0 otherwise
        X_g_tv=0;               %[boolean]    (1x1) 1 if X_g is time variant 0 otherwise
        X_tv=0;                 %[boolean]    (1x1) 1 if all the space-time varying coefficients are time variant 0 is they are all not
    end
    
    methods
        
        function obj = stem_data(stem_varset_g,stem_gridlist_g,stem_varset_r,stem_gridlist_r,stem_datestamp,shape,can_reset,stem_crossval,remote_correlated)
            %DESCRIPTION: is the constructor of the class stem_data
            %
            %INPUT
            %
            %stem_varset_g      - [stem_varset object]    (1x1) stem_varset object for the ground level variables
            %stem_gridlist_g    - [stem_gridlist object]  (1x1) stem_gridlist object for the ground level variables  
            %<stem_varset_r>    - [stem_varset object]    (1x1) (default: []) stem_varset object for the remote sensing variables
            %<stem_gridlist_r>  - [stem_gridlist object]  (1x1) (default: []) stem_gridlist object for the remote sensing variables
            %<stem_datestamp>   - [stem_datestamp object] (1x1) (default: []) stem_datestamp object with information on time steps
            %<shape>            - [struct]                      (default: world boundaries) geographic data structure loaded from a shapefile with the boundary of the geographic region
            %<can_reset>        - [boolean]               (1x1) (default: 0) if 1 the data are saved on disk and they can be reloaded using the method reset of this class after, for example, data transformation
            %<stem_crossval>    - [stem_crossval object]  (1x1) (default: []) stem_crossval object with information on crossvalidation
            %
            %OUTPUT
            %obj                - [stem_data object]      (1x1) the stem_data object
            if nargin<2
                error('Not enough input parameters');
            end
            obj.stem_varset_g=stem_varset_g;
            obj.stem_gridlist_g=stem_gridlist_g;
            if nargin==3
                error('stem_gridlist_r must be provided');
            end
            if nargin>2
                if not(isempty(stem_varset_r))
                    obj.stem_varset_r=stem_varset_r;
                    obj.stem_gridlist_r=stem_gridlist_r;
                end
            end
            if nargin>4
                obj.stem_datestamp=stem_datestamp;
            end

            if nargin>=6
                if not(isempty(shape))
                    obj.shape=shape;
                    obj.shape(1,1).Geometry='Line';
                else
                    obj.shape=shaperead('landareas.shp');
                end
            else
                obj.shape=shaperead('landareas.shp');
            end
            
            if nargin>=7
                if not(isempty(can_reset))
                    obj.can_reset=can_reset;
                else
                    obj.can_reset=0;
                end
            else
                obj.can_reset=0;
            end
            
            if nargin>=8
                if not(isempty(stem_crossval))
                    obj.stem_crossval=stem_crossval;
                    idx_var = obj.stem_varset_g.get_Y_index(obj.stem_crossval.variable_name);
                    if isempty(idx_var)
                        error('Cross-validation variable not found');
                    end
                end
            end
            
            if nargin>=9
                if not(isempty(remote_correlated))
                    obj.remote_correlated=remote_correlated;
                end
            end
            
%             if obj.can_reset
%                 obj.set_original;
%             end
            
            obj.update_data();
            obj.update_distance();
            if not(isempty(stem_varset_r))
                obj.update_M();
            end
            obj.remove_duplicated_sites();
        end
        
        function update_data(obj)
            %DESCRIPTION: generates the matrices Y,X_rg,X_beta,X_time,X_g 
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            disp('Generating data matrices...');
            %Y
            Y=[];
            for i=1:length(obj.stem_varset_g.Y)
                Y=cat(1,Y,obj.stem_varset_g.Y{i});
            end
            if not(isempty(obj.stem_varset_r))
                for i=1:length(obj.stem_varset_r.Y)
                    Y=cat(1,Y,obj.stem_varset_r.Y{i});
                end
            end
            obj.Y=Y;
            clear Y;
            %X_rg
            X_rg=[];
            if not(isempty(obj.stem_varset_r))
                if not(isempty(obj.stem_varset_r.X_rg))&&not(isempty(obj.stem_varset_g.X_rg))
                    done=0;
                    if size(obj.stem_varset_g.X_rg{1},3)==1 && size(obj.stem_varset_r.X_rg{1},3)==obj.T
                        for i=1:length(obj.stem_varset_g.X_rg)
                            X_rg=cat(1,X_rg,repmat(obj.stem_varset_g.X_rg{i},[1,1,obj.T]));
                        end
                        for i=1:length(obj.stem_varset_r.X_rg)
                            X_rg=cat(1,X_rg,obj.stem_varset_r.X_rg{i});
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_rg{1},3)==obj.T && size(obj.stem_varset_r.X_rg{1},3)==1 && not(done)
                        for i=1:length(obj.stem_varset_g.X_rg)
                            X_rg=cat(1,X_rg,obj.stem_varset_g.X_rg{i});
                        end
                        for i=1:length(obj.stem_varset_r.X_rg)
                            X_rg=cat(1,X_rg,repmat(obj.stem_varset_r.X_rg{i},[1,1,obj.T]));
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_rg{1},3)==size(obj.stem_varset_r.X_rg{1},3) && not(done)
                        for i=1:length(obj.stem_varset_g.X_rg)
                            X_rg=cat(1,X_rg,obj.stem_varset_g.X_rg{i});
                        end
                        for i=1:length(obj.stem_varset_r.X_rg)
                            X_rg=cat(1,X_rg,obj.stem_varset_r.X_rg{i});
                        end
                        done=1;
                    end
                end
                obj.X_rg=X_rg;
                if size(obj.X_rg,3)>1
                    obj.X_rg_tv=1;
                end
                clear X_rg;
            end
%             else
%                 for i=1:length(obj.stem_varset_g.X_rg)
%                     X_rg=cat(1,X_rg,obj.stem_varset_g.X_rg{i});
%                 end
%             end

            %X_beta
            X_beta=[];
            if not(isempty(obj.stem_varset_r))
                if not(isempty(obj.stem_varset_r.X_beta))&&not(isempty(obj.stem_varset_g.X_beta))
                    done=0;
                    if size(obj.stem_varset_g.X_beta{1},3)==1 && size(obj.stem_varset_r.X_beta{1},3)==obj.T
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i});
                            end
                            for i=1:length(obj.stem_varset_r.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_beta{1},3)==obj.T && size(obj.stem_varset_r.X_beta{1},3)==1 && not(done)
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_r.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_beta{i});
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_g.X_beta{1},3)==obj.T)&&(size(obj.stem_varset_r.X_beta{1},3)==obj.T) && not(done)
                        for t=1:size(obj.stem_varset_g.X_beta{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_r.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_g.X_beta{1},3)==1)&&(size(obj.stem_varset_r.X_beta{1},3)==1) && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_g.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i});
                        end
                        for i=1:length(obj.stem_varset_r.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_r.X_beta{i});
                        end
                        X_beta=X_temp;
                        done=1;
                    end
                else
                    if not(isempty(obj.stem_varset_g.X_beta))
                        done=0;
                        if size(obj.stem_varset_g.X_beta{1},3)==obj.T
                            for t=1:size(obj.stem_varset_g.X_beta{1},3)
                                X_temp=[];
                                for i=1:length(obj.stem_varset_g.X_beta)
                                    X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i}(:,:,t));
                                end
                                %X_temp=cat(1,X_temp,zeros(obj.stem_varset_r.N,size(X_temp,2),size(X_temp,3)));
                                X_beta(:,:,t)=X_temp;
                            end
                            done=1;
                        end
                        if size(obj.stem_varset_g.X_beta{1},3)==1 && not(done)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i});
                            end
                            %X_temp=cat(1,X_temp,zeros(obj.stem_varset_r.N,size(X_temp,2),size(X_temp,3)));
                            X_beta=X_temp;
                            done=1;
                        end
                    end
                end
            else
                if not(isempty(obj.stem_varset_g.X_beta))
                    done=0;
                    if size(obj.stem_varset_g.X_beta{1},3)==obj.T
                        for t=1:size(obj.stem_varset_g.X_beta{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_beta)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i}(:,:,t));
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_beta{1},3)==1 && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_g.X_beta)
                            X_temp=blkdiag(X_temp,obj.stem_varset_g.X_beta{i});
                        end
                        X_beta=X_temp;
                        done=1;
                    end
                end
            end
            obj.X_beta=X_beta;
            if size(obj.X_beta,3)>1
                obj.X_beta_tv=1;
            end
            clear X_beta;

            %X_time
            X_time=[];
            if not(isempty(obj.stem_varset_r))
                if not(isempty(obj.stem_varset_r.X_time))&&not(isempty(obj.stem_varset_g.X_time))
                    done=0;
                    if size(obj.stem_varset_g.X_time{1},3)==1 && size(obj.stem_varset_r.X_time{1},3)==obj.T
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i});
                            end
                            for i=1:length(obj.stem_varset_r.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_time{i}(:,:,t));
                            end
                            X_time(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_time{1},3)==obj.T && size(obj.stem_varset_r.X_time{1},3)==1 && not(done)
                        for t=1:obj.T
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_r.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_time{i});
                            end
                            X_time(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_g.X_time{1},3)==obj.T)&&(size(obj.stem_varset_r.X_time{1},3)==obj.T) && not(done)
                        for t=1:size(obj.stem_varset_g.X_time{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i}(:,:,t));
                            end
                            for i=1:length(obj.stem_varset_r.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_r.X_time{i}(:,:,t));
                            end
                            X_time(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if (size(obj.stem_varset_g.X_time{1},3)==1)&&(size(obj.stem_varset_r.X_time{1},3)==1) && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_g.X_time)
                            X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i});
                        end
                        for i=1:length(obj.stem_varset_r.X_time)
                            X_temp=blkdiag(X_temp,obj.stem_varset_r.X_time{i});
                        end
                        X_time=X_temp;
                        done=1;
                    end
                else
                    if not(isempty(obj.stem_varset_g.X_time))
                        done=0;
                        if size(obj.stem_varset_g.X_time{1},3)==obj.T
                            for t=1:size(obj.stem_varset_g.X_time{1},3)
                                X_temp=[];
                                for i=1:length(obj.stem_varset_g.X_time)
                                    X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i}(:,:,t));
                                end
                                X_time(:,:,t)=X_temp;
                            end
                            done=1;
                        end
                        if size(obj.stem_varset_g.X_time{1},3)==1 && not(done)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i});
                            end
                            X_time=X_temp;
                            done=1;
                        end
                        %X_time=cat(1,X_time,zeros(size(obj.Y,1)-size(X_time,1),size(X_time,2),size(X_time,3)));
                    end
                end
            else
                if not(isempty(obj.stem_varset_g.X_time))
                    done=0;
                    if size(obj.stem_varset_g.X_time{1},3)==obj.T
                        for t=1:size(obj.stem_varset_g.X_time{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_g.X_time)
                                X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i}(:,:,t));
                            end
                            X_time(:,:,t)=X_temp;
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_g.X_time{1},3)==1 && not(done)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_g.X_time)
                            X_temp=blkdiag(X_temp,obj.stem_varset_g.X_time{i});
                        end
                        X_time=X_temp;
                        done=1;
                    end
                end
            end
            if not(isempty(X_time))
                obj.X_time=X_time;
                if size(obj.X_time,3)>1
                    obj.X_time_tv=1;
                end
            end
            clear X_time

            %X_g
            if not(isempty(obj.stem_varset_g.X_g))
                X_g=[];
                for i=1:length(obj.stem_varset_g.X_g)
                    X_g=cat(1,X_g,obj.stem_varset_g.X_g{i});
                end
                obj.X_g=X_g;
                if size(obj.X_g,3)>1
                    obj.X_g_tv=1;
                end
                clear X_g;
            end
            obj.X_tv=obj.X_rg_tv | obj.X_time_tv | obj.X_g_tv;
            disp('Generation ended.');
        end
        
        function update_M(obj)
            %DESCRIPTION: generates the vector M
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT         
            %
            %none: the vector M is generated ad updated
            disp('Generating M replication vector...');
            %M
            M=[];
            blocks=[0 cumsum(obj.stem_varset_r.dim)];
            for j=1:obj.stem_varset_g.nvar
                dmax=distdim(distance(0,0,obj.stem_gridlist_r.grid{j}.pixel_side_w,obj.stem_gridlist_r.grid{j}.pixel_side_w), obj.stem_gridlist_g.grid{1}.unit, 'km');
                for i=1:size(obj.stem_gridlist_g.grid{j}.coordinate,1)
                    d=distdim(distance(obj.stem_gridlist_g.grid{j}.coordinate(i,:),obj.stem_gridlist_r.grid{j}.coordinate), obj.stem_gridlist_g.grid{1}.unit, 'km');
                    [m,idx]=min(d);
                    if d>dmax
                        warning(['Point ',num2str(i),' of ground level variable ',num2str(j),' does not belong to any pixel. The nearest pixel at ',num2str(m),' km is considered']);
                    end
                    M=[M;idx+blocks(j)];
                end
            end
            obj.M=M;
            disp('Generation ended.');
        end
        
        function update_distance(obj,type,force)
            %DESCRIPTION: generates the distance matrices
            %
            %INPUT
            %obj  - [stem_data object] (1x1) the stem_data object
            %<type> - [string]         (1x1) (Default: 'both') if type='ground' only the distance matrix for the ground level data is evaluated. 
                                                              %if type='remote' only the distance matrix for the remote sensing data is evaluated.
                                                              %if type='both' both the matrices are evaluated.
            %OUTPUT         
            %
            %none: the distance matrices are generated and updated
            
            %Distance matrices
            if nargin<2
                type='both';
            end
            if nargin<3
                force=0;
            end
            cmp=strcmp(type,{'both','ground','remote'});
            if sum(cmp)==0
                error('type must be ground, remote or both');
            end
            
            if strcmp(type,'ground')||strcmp(type,'both')
                if not(isempty(obj.stem_varset_g.X_g))||force
                    disp('Generating ground level distance matrices...');
                    obj.DistMat_g=obj.stem_gridlist_g.get_distance_matrix();
                    disp('Generation ended.');
                end
            end
            if strcmp(type,'remote')||strcmp(type,'both')
                if not(isempty(obj.stem_gridlist_r))&&not(isempty(obj.stem_varset_r.X_rg))
                    disp('Generating remote data distance matrices...');
                    obj.DistMat_r=obj.stem_gridlist_r.get_distance_matrix(obj.remote_correlated);
                    disp('Generation ended.');
                end
            end
        end
        
%         function google_map(obj,name,type)
%             if sum(strcmp(type,{'ground','remote'}))==0
%                 error('type must be either ground or remote');
%             end
%             if strcmp(type,'ground')
%                 %not supported yet
%             else
%                 if isempty(obj.stem_varset_r)
%                     disp('No remote variables in this model');
%                 else
%                     idx=find(strcmp(name,obj.stem_varset_r.Y_name));
%                     T=size(obj.stem_varset_r.Y{idx},2);
%                     datafile=[1,T,obj.stem_gridlist_r.grid{idx}.pixel_side_w];
%                     csvwrite('..\Data\google_bridge\parameters_kriging.csv',datafile);
%                     for t=1:T
%                         t
%                         min_value=nanmin(nanmin(obj.stem_varset_r.Y{idx}(:,t)));
%                         max_value=nanmax(nanmax(obj.stem_varset_r.Y{idx}(:,t)));
%                         data=obj.stem_varset_r.Y{idx}(:,t);
%                         value=(data-min_value)/(max_value-min_value);
%                         %color=stem_krig_result.toColor(value);
%                         datafile=[obj.stem_gridlist_r.grid{idx}.coordinate(:,1),obj.stem_gridlist_r.grid{idx}.coordinate(:,2),value,data];
%                         csvwrite(['..\Data\google_bridge\kriging',num2str(t),'.csv'],datafile);
%                     end
%                     winopen('..\Data\google_bridge\open_kriging.bat');
%                 end
%             end
%         end
        
        %data transform
        
        function detrend(obj)
            disp('Ground level data detrend started...');
            obj.stem_varset_g.detrend;
            disp('Ground level data detrend ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Remote sensing data detrend started...');
                obj.stem_varset_r.detrend;
                disp('Remote sensing data detrend ended.');
            end
            disp('Updtaing data matrices after detrend...');
            obj.update_data;            
        end
        
        function standardize_sbs(obj)
            disp('Ground level data site by site standardization started...');
            obj.stem_varset_g.standardize_sbs;
            disp('Ground level data site by site standardization ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Remote sensing data site by site standardization started...');
                obj.stem_varset_r.standardize_sbs;
                disp('Remote sensing data site by site standardization ended.');
            end
            disp('Updtaing data matrices after site by site standardization...');
            obj.update_data;
        end
        
        function standardize(obj)
            disp('Ground level data standardization started...');
            obj.stem_varset_g.standardize;
            disp('Ground level data standardization ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Remote sensing data standardization started...');
                obj.stem_varset_r.standardize;
                disp('Remote sensing data standardization ended.');
            end
            disp('Updtaing data matrices after standardization...');
            obj.update_data;
        end
        
        function log_transform(obj)
            disp('Ground level data log-transformation started...');
            obj.stem_varset_g.log_transform;
            disp('Ground level data log-transformation ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Remote sensing data log-transformation started...');
                obj.stem_varset_r.log_transform;
                disp('Remote sensing data log-transformation ended.');
            end
            disp('Updtaing data matrices after log-transformation...');
            obj.update_data;
        end
        
        function time_average(obj,n_steps)
            if nargin<2
                error('The subsampling factor must be provided');
            end
            if n_steps<=1
                error('n_steps must be greater than 1');
            end
            if round(n_steps)~=n_steps
                error('n_steps must be an integer value');
            end
   
            indices=0:n_steps:obj.T;
            if indices(end)~=obj.T
                indices=[indices,obj.T];
            end
            obj.stem_datestamp.average_stamps(indices);
            
            disp('Time averaging started...');
            for i=1:length(obj.stem_varset_g.Y)
                for j=1:length(indices)-1
                    Y_temp{i}(:,j)=nanmean(obj.stem_varset_g.Y{i}(:,indices(j)+1:indices(j+1)),2);
                end
            end
            obj.stem_varset_g.Y=Y_temp;
            clear Y_temp

            if not(isempty(obj.stem_varset_g.X_rg))
                if obj.stem_varset_g.X_rg_tv
                    for i=1:length(obj.stem_varset_g.X_rg)
                        for j=1:length(indices)-1
                            X_rg_temp{i}(:,:,j)=nanmean(obj.stem_varset_g.X_rg{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.stem_varset_g.X_rg=X_rg_temp;
                    clear X_rg_temp
                end
            end
            
            if not(isempty(obj.stem_varset_g.X_beta))
                if obj.stem_varset_g.X_beta_tv
                    for i=1:length(obj.stem_varset_g.X_beta)
                        for j=1:length(indices)-1
                            X_beta_temp{i}(:,:,j)=nanmean(obj.stem_varset_g.X_beta{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.stem_varset_g.X_beta=X_beta_temp;
                    clear X_beta_temp
                end
            end
            
%             %da attivare per stem
%             if not(isempty(obj.stem_varset_g.X_time))
%                 if obj.stem_varset_g.X_time_tv
%                     for i=1:length(obj.stem_varset_g.X_time)
%                         for j=1:length(indices)-1
%                             X_time_temp{i}(:,:,j)=nanmean(obj.stem_varset_g.X_time{i}(:,:,indices(j)+1:indices(j+1)),3);
%                         end
%                     end
%                     obj.stem_varset_g.X_time=X_time_temp;
%                     clear X_time_temp
%                 end
%             end   
            
            if not(isempty(obj.stem_varset_g.X_g))
                if obj.stem_varset_g.X_g_tv
                    for i=1:length(obj.stem_varset_g.X_g)
                        for j=1:length(indices)-1
                            X_g_temp{i}(:,:,j,:)=nanmean(obj.stem_varset_g.X_g{i}(:,:,indices(j)+1:indices(j+1),:),3);
                        end
                    end
                    obj.stem_varset_g.X_g=X_g_temp;        
                    clear X_g_temp
                end
            end
            
            
            if not(isempty(obj.stem_varset_r))
                for i=1:length(obj.stem_varset_r.Y)
                    for j=1:length(indices)-1
                        Y_temp{i}(:,j)=nanmean(obj.stem_varset_r.Y{i}(:,indices(j)+1:indices(j+1)),2);
                    end
                end
                obj.stem_varset_r.Y=Y_temp;
                clear Y_temp
                
                if not(isempty(obj.stem_varset_g.X_rg))
                    if obj.stem_varset_r.X_rg_tv
                        for i=1:length(obj.stem_varset_r.X_rg)
                            for j=1:length(indices)-1
                                X_rg_temp{i}(:,:,j)=nanmean(obj.stem_varset_r.X_rg{i}(:,:,indices(j)+1:indices(j+1)),3);
                            end
                        end
                        obj.stem_varset_r.X_rg=X_rg_temp;
                        clear X_rg_temp
                    end
                end
                
                if not(isempty(obj.stem_varset_r.X_beta))
                    if obj.stem_varset_r.X_beta_tv
                        for i=1:length(obj.stem_varset_r.X_beta)
                            for j=1:length(indices)-1
                                X_beta_temp{i}(:,:,j)=nanmean(obj.stem_varset_r.X_beta{i}(:,:,indices(j)+1:indices(j+1)),3);
                            end
                        end
                        obj.stem_varset_r.X_beta=X_beta_temp;
                        clear X_beta_temp
                    end
                end
                
%                 %da attivare per stem
%                 if not(isempty(obj.stem_varset_r.X_time))
%                     if obj.stem_varset_r.X_time_tv
%                         for i=1:length(obj.stem_varset_r.X_time)
%                             for j=1:length(indices)-1
%                                 X_time_temp{i}(:,:,j)=nanmean(obj.stem_varset_r.X_time{i}(:,:,indices(j)+1:indices(j+1)),3);
%                             end
%                         end
%                         obj.stem_varset_r.X_time=X_time_temp;
%                     end
%                 end
                
            end

            obj.update_data;    
            disp('Time averaging ended.');
        end
        
        function time_crop(obj,dates_or_indices)
            % cut data in time
            % dates_or_indices: can be either a {2x1} cell vector of starting and ending date in the format dd-mm-yyyy or a T'x1 vector of (possibly non consecutive) indices. T'<=T.
            if iscell(dates_or_indices)
                start_date=datenum(dates_or_indices{1},'dd-mm-yyyy');
                end_date=datenum(dates_or_indices{2},'dd-mm-yyyy');
                if start_date>end_date
                    error('Starting date must be lower or equal to ending date');
                end
                sd=find(obj.stem_datestamp.stamp==start_date);
                ed=find(obj.stem_datestamp.stamp==end_date);
                if isempty(sd)||isempty(ed)
                    error('Dates out of range');
                end
                indices=sd:ed;
            else
                if min(dates_or_indices)<0
                    error('Starting date must be higher than zero');
                end
                if max(dates_or_indices)>obj.T
                    error('Ending date must be lower or equal to the total number of temporal steps');
                end
                indices=dates_or_indices;
            end
            
            disp('Time crop started...');
            Y=obj.stem_varset_g.Y;
            for i=1:length(Y)
                Y{i}=Y{i}(:,indices);
            end
            obj.stem_varset_g.Y=Y;
            if not(isempty(obj.stem_varset_g.X_beta))
                if obj.stem_varset_g.X_beta_tv
                    X_beta=obj.stem_varset_g.X_beta;
                    for i=1:length(X_beta)
                        X_beta{i}=X_beta{i}(:,:,indices);
                    end
                    obj.stem_varset_g.X_beta=X_beta;
                end
            end
            if not(isempty(obj.stem_varset_g.X_rg))
                if obj.stem_varset_g.X_rg_tv
                    X_rg=obj.stem_varset_g.X_rg;
                    for i=1:length(X_rg)
                        X_rg{i}=X_rg{i}(:,:,indices);
                    end
                    obj.stem_varset_g.X_rg=X_rg;
                end
            end               
            if not(isempty(obj.stem_varset_g.X_time))
                if obj.stem_varset_g.X_time_tv
                    X_time=obj.stem_varset_g.X_time;
                    for i=1:length(X_time)
                        X_time{i}=X_time{i}(:,:,indices);
                    end
                    obj.stem_varset_g.X_time=X_time;
                end
            end   
            if not(isempty(obj.stem_varset_g.X_g))
                if obj.stem_varset_g.X_g_tv
                    X_g=obj.stem_varset_g.X_g;
                    for i=1:length(X_g)
                        X_g{i}=X_g{i}(:,:,indices,:);
                    end
                    obj.stem_varset_g.X_g=X_g;
                end
            end
            
            if not(isempty(obj.stem_varset_r))
                if not(isempty(obj.stem_varset_r.Y))
                    Y=obj.stem_varset_r.Y;
                    for i=1:length(Y)
                        Y{i}=Y{i}(:,indices);
                    end
                    obj.stem_varset_r.Y=Y;
                    if not(isempty(obj.stem_varset_r.X_beta))
                        if obj.stem_varset_r.X_beta_tv
                            X_beta=obj.stem_varset_r.X_beta;
                            for i=1:length(obj.stem_varset_r.X_beta)
                                X_beta{i}=X_beta{i}(:,:,indices);
                            end
                            obj.stem_varset_r.X_beta=X_beta;
                        end
                    end
                    if not(isempty(obj.stem_varset_r.X_rg))
                        if obj.stem_varset_r.X_rg_tv
                            X_rg=obj.stem_varset_r.X_rg;
                            for i=1:length(X_rg)
                                X_rg{i}=X_rg{i}(:,:,indices);
                            end
                            obj.stem_varset_r.X_time=X_rg;
                        end
                    end                    
                    if not(isempty(obj.stem_varset_r.X_time))
                        if obj.stem_varset_r.X_time_tv
                            X_time=obj.stem_varset_r.X_time;
                            for i=1:length(X_time)
                                X_time{i}=X_time{i}(:,:,indices);
                            end
                            obj.stem_varset_r.X_time=X_time;
                        end
                    end
                end
            end
            
            obj.stem_datestamp.subset_stamps(indices);
            %looks for line of all missing for the sparse grids of the ground level data
            changed=0;
            for i=1:obj.stem_varset_g.nvar
                indices=sum(isnan(obj.stem_varset_g.Y{i}),2)==size(obj.stem_varset_g.Y{i},2);
                if sum(indices)>0
                    obj.stem_varset_g.Y{i}(indices,:)=[];
                    if not(isempty(obj.stem_varset_g.X_beta))
                        obj.stem_varset_g.X_beta{i}(indices,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_g.X_time))
                        obj.stem_varset_g.X_time{i}(indices,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_g.X_g))
                        obj.stem_varset_g.X_g{i}(indices,:,:,:)=[];
                    end
                    if not(isempty(obj.stem_varset_g.X_rg))
                        obj.stem_varset_g.X_rg{i}(indices,:,:)=[];
                    end
                    obj.stem_gridlist_g.grid{i}.coordinate(indices,:)=[];
                    disp(['Deleted ',num2str(sum(indices)),' site(s) for the ground variable ',obj.stem_varset_g.Y_name{i},' due to all missing.']);
                    changed=1;
                end
            end
            if changed
                disp('Updating ground distance matrix after time crop...');
                obj.update_distance('ground'); %only ground because the remote data are not deleted from the data matrix even if they are NaN for all th time steps
                disp('Update ended.');
                if not(isempty(obj.stem_varset_r))
                    disp('Updating M replication vector after time crop...');
                    obj.update_M;
                    disp('Update ended.');
                end
            end
            disp('Updating data matrix after time crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Time crop ended.');
        end     
        
        function space_crop(obj,box)
            % cut data in space
            lat_min=box(1);
            lat_max=box(2);
            lon_min=box(3);
            lon_max=box(4);
            if nargin<2
                error('The bounding-box must be provided');
            end
            if lat_min>lat_max
                error('The lat_min value must be lower than the lat_max value');
            end
            if lon_min>lon_max
                error('The lon_min value must be lower than the lon_max value');
            end
            disp('Ground level data space crop started...');
            for i=1:obj.stem_varset_g.nvar
                GY=obj.stem_gridlist_g.grid{i}.coordinate;
                indices = GY(:,2) >= lon_min & ...
                    GY(:,2) <= lon_max & ...
                    GY(:,1) >= lat_min & ...
                    GY(:,1) <= lat_max;
                if sum(indices)>0
                    if strcmp(obj.stem_gridlist_g.grid{i}.grid_type,'regular')
                        grid_lat=obj.stem_gridlist_g.grid{i}.coordinate(:,1);
                        grid_lon=obj.stem_gridlist_g.grid{i}.coordinate(:,2);
                        grid_lat=reshape(grid_lat,obj.stem_gridlist_g.grid{i}.grid_size);
                        grid_lat=grid_lat(:,1);
                        grid_lon=reshape(grid_lon,obj.stem_gridlist_g.grid{i}.grid_size);
                        grid_lon=grid_lon(1,:);
                        grid_lat(grid_lat<lat_min)=[];
                        grid_lat(grid_lat>lat_max)=[];
                        grid_lon(grid_lon<lon_min)=[];
                        grid_lon(grid_lon>lon_max)=[];
                    end
                    obj.stem_gridlist_g.grid{i}.coordinate=obj.stem_gridlist_g.grid{i}.coordinate(indices,:);
                    
                    if strcmp(obj.stem_gridlist_g.grid{i}.grid_type,'regular')
                        obj.stem_gridlist_g.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                    end                    

                    obj.stem_varset_g.Y{i}=obj.stem_varset_g.Y{i}(indices,:);
                    
                    if not(isempty(obj.stem_varset_g.X_rg))
                        obj.stem_varset_g.X_rg{i}=obj.stem_varset_g.X_rg{i}(indices,:,:);
                    end
                    
                    if not(isempty(obj.stem_varset_g.X_beta))
                        obj.stem_varset_g.X_beta{i}=obj.stem_varset_g.X_beta{i}(indices,:,:);
                    end
                    if not(isempty(obj.stem_varset_g.X_time))
                        obj.stem_varset_g.X_time{i}=obj.stem_varset_g.X_time{i}(indices,:,:);
                    end
                    if not(isempty(obj.stem_varset_g.X_g))
                        obj.stem_varset_g.X_g{i}=obj.stem_varset_g.X_g{i}(indices,:,:,:);
                    end
                else
                    error(['    Variable ',obj.stem_varset_g.Y_name{i},' does not have sites in the crop area.']);
                end
            end
            disp('Ground level data space crop ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Remote sensing data space crop started...');
                for i=1:obj.stem_varset_r.nvar
                    GY=obj.stem_gridlist_r.grid{i}.coordinate;
                    indices = GY(:,2) >= lon_min & ...
                        GY(:,2) <= lon_max & ...
                        GY(:,1) >= lat_min & ...
                        GY(:,1) <= lat_max;
                    
                    if sum(indices)>0
                        if strcmp(obj.stem_gridlist_r.grid{i}.grid_type,'regular')
                            grid_lat=obj.stem_gridlist_r.grid{i}.coordinate(:,1);
                            grid_lon=obj.stem_gridlist_r.grid{i}.coordinate(:,2);
                            grid_lat=reshape(grid_lat,obj.stem_gridlist_r.grid{i}.grid_size);
                            grid_lat=grid_lat(:,1);
                            grid_lon=reshape(grid_lon,obj.stem_gridlist_r.grid{i}.grid_size);
                            grid_lon=grid_lon(1,:);
                            grid_lat(grid_lat<lat_min)=[];
                            grid_lat(grid_lat>lat_max)=[];
                            grid_lon(grid_lon<lon_min)=[];
                            grid_lon(grid_lon>lon_max)=[];
                        end
                       
                        obj.stem_gridlist_r.grid{i}.coordinate=obj.stem_gridlist_r.grid{i}.coordinate(indices,:);
                        if strcmp(obj.stem_gridlist_r.grid{i}.grid_type,'regular')
                            obj.stem_gridlist_r.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                        end
                        obj.stem_varset_r.Y{i}=obj.stem_varset_r.Y{i}(indices,:);
                        if not(isempty(obj.stem_varset_r.X_rg))
                            obj.stem_varset_r.X_rg{i}=obj.stem_varset_r.X_rg{i}(indices,:,:);
                        end
                        if not(isempty(obj.stem_varset_r.X_beta))
                            obj.stem_varset_r.X_beta{i}=obj.stem_varset_r.X_beta{i}(indices,:,:);
                        end
                        if not(isempty(obj.stem_varset_r.X_time))
                            obj.stem_varset_r.X_time{i}=obj.stem_varset_r.X_time{i}(indices,:,:);
                        end
                    else
                        error(['    Variable ',obj.stem_varset_r.Y_name{i},' does not have sites in the crop area.']);
                    end
                end
                disp('Remote sensing data space crop ended.');
            end
            disp('Updating data matrices after space crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after space crop...');
            obj.update_distance;
            disp('Update ended.');

            if not(isempty(obj.stem_varset_r))
                disp('Updating M replication vector after space crop...');
                obj.update_M;
                disp('Update ended.');
            end

        end   
        
        function site_crop(obj,type,var_name,indices)
            if sum(strcmp(type,{'ground','remote'}))==0
                error('Type must be either ground or remote');
            end
            if min(indices<1)
                error('The minimum value of indices cannot be lower than 1');
            end
            if strcmp(type,'ground')
                idx_var=obj.stem_varset_g.get_Y_index(var_name);
                if isempty(idx_var)
                    error('Variable not found');
                end
                N=obj.stem_varset_g.N;
                if max(indices>N)
                    error(['The maximum value of indices cannot be greater than ',num2str(N)]);
                end
                obj.stem_gridlist_g.grid{idx_var}.coordinate(indices,:)=[];
                obj.stem_varset_g.Y{idx_var}(indices,:)=[];
                
                if not(isempty(obj.stem_varset_g.X_rg))
                    obj.stem_varset_g.X_rg{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_g.X_beta))
                    obj.stem_varset_g.X_beta{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_g.X_time))
                    obj.stem_varset_g.X_time{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_g.X_g))
                    obj.stem_varset_g.X_g{idx_var}(indices,:,:,:)=[];
                end
            else
                if isempty(obj.stem_varset_r)
                    error('No remote data');
                end
                idx_var=obj.stem_varset_r.get_Y_index(var_name);
                if isempty(idx_var)
                    error('Variable not found');
                end
                obj.stem_gridlist_r.grid{idx_var}.coordinate(indices,:)=[];
                obj.stem_varset_r.Y{idx_var}(indices,:)=[];
                if not(isempty(obj.stem_varset_r.X_beta))
                    obj.stem_varset_r.X_beta{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.stem_varset_r.X_time))
                    obj.stem_varset_r.X_time{idx_var}(indices,:,:)=[];
                end
            end
            disp('Updating data matrices after site crop...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after site crop...');
            obj.update_distance(type);
            disp('Update ended.');
            if not(isempty(obj.stem_varset_r))
                disp('Updating M replication vector after site crop...');
                obj.update_M;
                disp('Update ended.');
            end
        end
        
        function remove_duplicated_sites(obj)
            %toll is the distance below which two sites are considered
            %colocated
            
            disp('Look for duplicated sites...')
            for i=1:obj.stem_varset_g.nvar
                idx=obj.stem_gridlist_g.grid{i}.duplicated_sites;
                if not(isempty(idx))
                    disp(['Removing ',num2str(length(idx)),' replicated sites for ground variable ',obj.stem_varset_g.Y_name{i}]);
                    obj.site_crop('ground',obj.stem_varset_g.Y_name{i},idx);
                    obj.stem_gridlist_g.grid{i}.duplicated_sites=[];
                end
            end
            
            if not(isempty(obj.stem_varset_r))
                for i=1:obj.stem_varset_r.nvar
                    idx=obj.stem_gridlist_r.grid{i}.duplicated_sites;
                    if not(isempty(idx))
                        disp(['Removing ',num2str(length(idx)),' replicated sites for remote variable ',obj.stem_varset_r.Y_name{i}]);
                        obj.site_crop('remote',obj.stem_varset_r.Y_name{i},idx);
                        obj.stem_gridlist_r.grid{i}.duplicated_sites=[];
                    end
                end
            end
            disp('Operation ended.')
        end        
        
        %export functions
        function N = N(obj)
            %return the value of N from the stem_varset_g and the stem_varset_r objects
            N=obj.stem_varset_g.N;
            if not(isempty(obj.stem_varset_r))
                N=N+obj.stem_varset_r.N;
            end
        end
        
        function Nr = Nr(obj)
            if not(isempty(obj.stem_varset_r))
                Nr = obj.stem_varset_r.N;
            else
                Nr=0;
            end
        end
        
        function Ng = Ng(obj)
            Ng = obj.stem_varset_g.N;
        end
        
        function T = T(obj)
            %return the value of T from the stem_varset_g and the stem_varset_r objects
            T=obj.stem_varset_g.T;
        end
        
        function nvar = nvar(obj)
            %return the value of nvar from the stem_varset_g and the stem_varset_r objects
            nvar=obj.stem_varset_g.nvar;
            if not(isempty(obj.stem_varset_r))
                nvar=nvar+obj.stem_varset_r.nvar;
            end
        end
        
        function dim = dim(obj)
            %return the value of dim from the stem_varset_g and the stem_varset_r objects
            dim=obj.stem_varset_g.dim;
            if not(isempty(obj.stem_varset_r))
                dim=[dim obj.stem_varset_r.dim];
            end
        end
        
        %set methods
        
        function set.stem_varset_g(obj,stem_varset_g)
           if not(isa(stem_varset_g,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           obj.stem_varset_g=stem_varset_g;
        end
        
        function set.stem_varset_r(obj,stem_varset_r)
           if not(isa(stem_varset_r,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           
           if not(length(stem_varset_r.dim)==length(obj.stem_varset_g.dim))
               error('stem_varset_r must contain the same number of variables of stem_varset_g');
           end
           
           if not(size(stem_varset_r.Y{1},2)==size(obj.stem_varset_g.Y{1},2))
               error('The number of temporal steps cannot differ between stem_varset_r and stem_varset_g');
           end
           
           if not(isempty(stem_varset_r.X_g))
               error('X_g must be empty in stem_varset_r');
           end
          
           obj.stem_varset_r=stem_varset_r;
        end        
        
        function set.stem_gridlist_g(obj,stem_gridlist_g)
            if not(isa(stem_gridlist_g,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            if not(length(stem_gridlist_g.grid)==length(obj.stem_varset_g.Y))
                error('The number of stem_grids must be equal to the q');
            end
            for i=1:length(stem_gridlist_g.grid)
                if not(size(stem_gridlist_g.grid{i}.coordinate,1)==size(obj.stem_varset_g.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(stem_gridlist_g.grid{i}.site_type,'point'))
                    error('Only point data are supported in stem_gridlist_g');
                end
            end
            obj.stem_gridlist_g=stem_gridlist_g;
        end
        
        function set.stem_gridlist_r(obj,stem_gridlist_r)
            if not(isa(stem_gridlist_r,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            if not(length(stem_gridlist_r.grid)==length(obj.stem_varset_r.Y))
                error('The number of stem_grids must be equal to the q');
            end            
            for i=1:length(stem_gridlist_r.grid)
                if not(size(stem_gridlist_r.grid{i}.coordinate,1)==size(obj.stem_varset_r.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(stem_gridlist_r.grid{i}.site_type,'pixel'))
                    error('The grids of stem_gridlist_r must be grids of pixels');
                end
                if not(strcmp(stem_gridlist_r.grid{i}.pixel_shape,'square'))
                    error('Only square pixels are supported. Check pixel shape');
                end
                if not(stem_gridlist_r.grid{i}.pixel_side_w==stem_gridlist_r.grid{i}.pixel_side_h)
                    error('Only square pixels are supported. Check pixel_side_w and pixel_side_h');
                end
            end 
            if not(strcmp(stem_gridlist_r.grid{1}.unit,obj.stem_gridlist_g.grid{1}.unit))
                error('Both the stem_gridlist objects must contain grids with the same unit');
            end
            obj.stem_gridlist_r=stem_gridlist_r;
        end  

        function set.stem_crossval(obj,stem_crossval)
            if not(isa(stem_crossval,'stem_crossval'))
                error('stem_crossval must be of class stem_crossval');
            end
            obj.stem_crossval=stem_crossval;
        end
        
        function set.stem_datestamp(obj,stem_datestamp)
            if not(isa(stem_datestamp,'stem_datestamp'))
                error('stem_datestamp must be a stem_datestamp object');
            end
            if not(obj.stem_varset_g.T==length(stem_datestamp.stamp))
                error('The number of datestamps differs from T');
            end
            obj.stem_datestamp=stem_datestamp;
        end

        function set.shape(obj,shape)
            obj.shape=shape;
        end
        

       
%         function remove_var(obj,indices)
%             indices=sort(indices);
%             if (min(indices)<1)||(max(indices)>obj.stem_varset.nvar)
%                 error('Indices out of bounds');
%             end
%             if length(indices)>1
%                 temp1=indices';
%                 temp2=circshift(temp1,1);
%                 temp=temp1-temp2;
%                 if nnz(temp)~=length(temp)
%                     error('Indices cannot be replicated');
%                 end
%             end
%             disp('Variable subsetting started...');
%             for i=1:length(indices)
%                 disp(['    Variable ',obj.stem_varset.name{indices(i)},' removed']);
%             end
%             obj.stem_varset.Y(indices)=[];
%             obj.stem_varset.name(indices)=[];
%             obj.Y=obj.stem_varset.get_Y;
%             if not(isempty(obj.stem_covset))
%                 obj.stem_covset.X(indices)=[];
%                 obj.stem_covset.name(indices)=[];
%                 obj.X=obj.stem_covset.get_X;
%             end
%             obj.stem_gridlist.grid(indices)=[];
%             obj.DistMat = obj.get_distance_matrix(obj.stem_gridlist); 
%             disp('Variable subsetting ended.');
%         end
        
%         function remove_cov(obj,variable_name,covariate_name)
%             if not(isempty(obj.stem_covset))
%                 if nargin<3
%                     error('All inputs must be provided');
%                 end
%                 if isempty(variable_name)
%                     index_var=1:obj.stem_varset.nvar;
%                 else
%                     index_var=obj.stem_varset.get_index(variable_name);
%                 end
%                 if not(isempty(index_var))
%                     for i=index_var
%                         index_cov=obj.stem_covset.get_index(i,covariate_name);
%                         if not(isempty(index_cov))
%                             obj.stem_covset.X{i}(:,index_cov,:)=[];
%                             obj.stem_covset.name{i}(index_cov)=[];
%                         else
%                             error('Covariate name not recognized for at least one variable');
%                         end
%                     end
%                 else
%                     error('Variable name not recognized');
%                 end
%                 obj.X=obj.stem_covset.get_X;
%             else
%                 warning('No covariates in stem_data object');
%             end
%         end
        
%         function reset(obj)
%             % restore original data
%             if obj.can_reset
%                 disp('Restoring original data...');
% 
%                 disp('Original data restored.');
%             else
%                 warning('Reset is not allow. Data must be reloaded by using the class constructor');
%             end
%         end
        
%         function print(obj)
%             % print data properties
%             disp('*********************************');
%             disp('*    STEM3 Data description     *');
%             disp('*********************************');
%             disp(' ');
%             disp(['Number of variables: ',num2str(obj.stem_varset.nvar)]);
%             disp(' ');
%             disp('Bounding box:');
%             disp(['  Latitude min : ',num2str(obj.stem_gridlist.box(1)),'°'])
%             disp(['  Latitude max : ',num2str(obj.stem_gridlist.box(2)),'°'])
%             disp(['  Longitude min: ',num2str(obj.stem_gridlist.box(3)),'°'])
%             disp(['  Longitude max: ',num2str(obj.stem_gridlist.box(4)),'°'])
%             disp(' ');
%             disp('Date and time:')
%             disp(['  Stating date  : ',datestr(obj.stem_datestamp.date_start)]);
%             disp(['  Ending date   : ',datestr(obj.stem_datestamp.date_end)]);
%             disp(['  Temporal steps: ',num2str(obj.stem_varset.T)]);
%             disp(' ');
%             disp('Variable description:');
%             output{1,1}='NAME';
%             output{1,2}='# COVAR.';
%             output{1,3}='# SITES';
%             output{1,4}='GRID SIZE';
%             output{1,5}='GRID TYPE';
%             output{1,6}='SITE TYPE';
%             output{1,7}='NaN';
%             output{1,8}='MEAN';
%             output{1,9}='STD';
%             for i=1:obj.stem_varset.nvar
%                 output{i+1,1}=obj.stem_varset.name{i};
%                 if not(isempty(obj.stem_covset))
%                     output{i+1,2}=num2str(obj.stem_covset.nbeta(i));
%                 else
%                     output{i+1,2}='0';
%                 end
%                 output{i+1,3}=num2str(obj.stem_varset.dim(i));
%                 if not(isempty(obj.stem_gridlist.grid{i}.grid_size))
%                     output{i+1,4}=[num2str(obj.stem_gridlist.grid{i}.grid_size(1)),'x',num2str(obj.stem_gridlist.grid{i}.grid_size(2))];
%                 else
%                     output{i+1,4}='/';    
%                 end
%                 output{i+1,5}=obj.stem_gridlist.grid{i}.grid_type;
%                 output{i+1,6}=obj.stem_gridlist.grid{i}.site_type;
%                 output{i+1,7}=[num2str(sum(isnan(obj.stem_varset.Y{i}(:)))/(obj.stem_varset.dim(i)*obj.stem_varset.T)*100,'%5.2f'),'%'];
%                 output{i+1,8}=num2str(nanmean(obj.stem_varset.Y{i}(:)),'%5.5f');
%                 output{i+1,9}=num2str(nanstd(obj.stem_varset.Y{i}(:)),'%5.5f');
%             end
%             disp(output)
%             disp(' ');
%             if not(isempty(obj.stem_covset))
%                 for i=1:obj.stem_varset.nvar
%                     output=[];
%                     if obj.stem_covset.nbeta(i)>0
%                         disp(['Covariate description for variable ',obj.stem_varset.name{i},':'])
%                         output{1,1}='NAME';
%                         output{1,2}='MEAN';
%                         output{1,3}='STD';
%                         for j=1:obj.stem_covset.nbeta(i)
%                             temp=obj.stem_covset.X{i}(:,j,:);
%                             output{j+1,1}=obj.stem_covset.name{i}{(j)};
%                             output{j+1,2}=num2str(mean(temp(:)));
%                             output{j+1,3}=num2str(std(temp(:)));
%                         end
%                     end
%                     disp(output);
%                 end
%             end
%             if obj.stem_varset.log_transformed
%                 disp(' ');
%                 disp('Variables are log-transformed. (Use the stem_model.data_reset() method to load the original data.)');
%             end
%             if obj.stem_varset.standardized
%                 disp(' ');
%                 disp('Variables are standardized. (Use the stem_model.data_reset() method to load the original data.)');
%             end
%             if not(isempty(obj.stem_covset))
%                 if obj.stem_covset.log_transformed
%                     disp(' ');
%                     disp('Covariates are log-transformed. (Use the stem_model.data_reset() method to load the original data.)');
%                 end
%                 if obj.stem_covset.standardized
%                     disp(' ');
%                     disp('Covariates are standardized. (Use the stem_model.data_reset() method to load the original data.)');
%                 end
%             end
%         end
        
%         function [v,theta,sigma_eps] = variogram(obj,n_steps,time_steps,graph)
%             if nargin<2
%                 n_steps=50;
%             end
%             if nargin<3
%                 time_steps=1:obj.stem_varset.T;
%             end
%             if nargin<4
%                 graph=0;
%             end
%             if isempty(time_steps)
%                 time_steps=1:obj.stem_varset.T;
%             end
%             options = optimset('Display','off');
%             
%             v=zeros(length(obj.stem_varset.dim));
%             theta=zeros(length(obj.stem_varset.dim));
%             sigma_eps=zeros(length(obj.stem_varset.dim));
%             
%             blocks=[0 cumsum(obj.stem_varset.dim)];
%             for k=1:obj.stem_varset.nvar
%                 for h=k:obj.stem_varset.nvar 
%                     y_sub1=obj.stem_varset.Y{h};
%                     y_sub2=obj.stem_varset.Y{k};
%                     DistMat_sub=obj.DistMat(blocks(h)+1:blocks(h+1),blocks(k)+1:blocks(k+1));
%                     max_distance=max(DistMat_sub(:));
%                     step=max_distance/n_steps;
%                     groups=0:step:max_distance;
%                     group=zeros(length(groups)-1,1);
%                     count=zeros(length(groups)-1,1);
%                     
%                     if h~=k
%                         tot=obj.stem_varset.dim(k)*obj.stem_varset.dim(h)*length(time_steps);
%                     else
%                         tot=((obj.stem_varset.dim(k)*(obj.stem_varset.dim(k)+1))/2)*length(time_steps);
%                     end
%                     vx=zeros(tot,1);
%                     vy=zeros(tot,1);
%                     index=1;
%                     for t=time_steps
%                         if h~=k
%                             y1=y_sub1(:,t);
%                             y2=y_sub2(:,t);
%                             for i=1:length(y1)
%                                 d=DistMat_sub(i,:);
%                                 sq=0.5*(y1(i)-y2).^2;
%                                 vx(index:index+length(d)-1)=d;
%                                 vy(index:index+length(d)-1)=sq;
%                                 index=index+length(d);
%                             end
%                         else
%                             y1=y_sub1(:,t);
%                             y2=y_sub2(:,t);
%                             for i=1:length(y1)
%                                 d=DistMat_sub(i,i:end);
%                                 sq=0.5*(y1(i)-y2(i:end)).^2;
%                                 vx(index:index+length(d)-1)=d;
%                                 vy(index:index+length(d)-1)=sq;
%                                 index=index+length(d);
%                             end
%                         end
%                     end
%                             
%                     for i=2:length(groups)
%                         data=vy(vx>=groups(i-1)&vx<groups(i));
%                         group(i-1)=nanmean(data);
%                         count(i-1)=length(data(isnotnan(data)));
%                     end
%                     xv=(groups(1:10)+(step/2))';
%                     yv=group(1:10);
%                     xv=xv(not(isnan(yv)));
%                     yv=yv(not(isnan(yv)));
%                     b=regress(yv,[ones(length(xv),1) xv]);
%                     v0=b(1);
%                     if v0<0
%                         v0=0.7;
%                     end
%                             
%                     L=isnotnan(group);
%                     L(1)=0;
%                     g=groups(1:end-1)+(step/2);
%                     x1=g(L)';
%                     y1=group(L);
%                     weights=count(L);
%                                                 
%                     if 0
%                         x0=[1 20];
%                         f = @(x,xdata) x(1)*(1-exp(-xdata/x(2)))+v0;'x';'xdata';
%                         c = lsqcurvefit(f,x0,x1,y1,[],[],options);
%                         theta0(index_m)=c(2);
%                         vv(index_m)=c(1);
%                         index_m=index_m+1;
%                     else
%                         o = fitoptions('Method','NonlinearLeastSquares');
%                         o.Weights=weights;
%                         o.StartPoint=[1 max_distance/10 v0];
%                         f = fittype('a*(1-exp(-x/b))+c');
%                         try
%                             %fit1 = fit(x1,y1-v0,f,o);
%                             fit1 = fit(x1,y1,f,o);
%                             c=coeffvalues(fit1);
%                         catch
%                             disp('Error');
%                         end
%                     end
%                             
%                     if graph
%                         figure
%                         xx=0:0.1:max(groups);
%                         yy=c(1)*(1-exp(-xx./c(2)))+c(3);
%                         plot(groups(1:end-1)+(step/2),group,'.');
%                         hold on
%                         plot(xx,yy,'-r');
%                         xlabel('Distance');
%                         if h==k
%                             title(['Variogram of ',obj.stem_varset.name{k}]);
%                         else
%                             title(['Cross-variogram between ',obj.stem_varset.name{k},' and ',obj.stem_varset.name{h}]);
%                         end
%                     end
%                 end
%             end
%         end
        
%         function plot_cov(obj,variable_name,covariate_name,trend_type)
%             % plot covariates against the specified variable
%             % variable_name: variable name string
%             if nargin<3
%                 covariate_name=[];
%             end
%             if nargin<4
%                 trend_type='linear';
%             end
%             vindex=obj.stem_varset.get_index(variable_name);
%             if not(isempty(vindex))
%                 if isempty(covariate_name)
%                     cov_indices=1:obj.stem_covset.nbeta(vindex);
%                     l=obj.stem_covset.nbeta(vindex)^0.5;
%                     if round(l)^2==obj.stem_covset.nbeta(vindex)
%                         rows=l;
%                         cols=l;
%                     else
%                         rows=ceil(l);
%                         cols=round(l);
%                     end
%                     temp=0;
%                 else
%                     cindex=obj.stem_covset.get_index(vindex,covariate_name);
%                     cov_indices=cindex;
%                     rows=1;
%                     cols=1;
%                     temp=cindex-1;
%                 end
% 
%                 for i=cov_indices
%                     subplot(rows,cols,i-temp);
%                     a=obj.stem_covset.X{vindex}(:,i,:);
%                     b=obj.stem_varset.Y{vindex};
%                     a=a(:);
%                     b=b(:);
%                     L=isnotnan(b);
%                     a=a(L);
%                     b=b(L);
%                     if strcmp(trend_type,'linear')
%                         coeff = regress(b,[ones(length(a),1),a]);
%                         xplot=min(a):range(a)/200:max(a);
%                         yplot=coeff(1)+xplot.*coeff(2);
%                     else
%                         coeff = regress(b,[ones(length(a),1),a,a.^2]);
%                         xplot=min(a):range(a)/200:max(a);
%                         yplot=coeff(1)+xplot.*coeff(2)+(xplot.^2).*coeff(3);
%                     end
%                     plot(a,b,'.b');
%                     hold on
%                     plot(xplot,yplot,'-k','LineWidth',2);
%                     ylabel(variable_name);
%                     xlabel(obj.stem_covset.name{vindex}{i});
%                 end
%             else
%                 error('Variable name not recognized');
%             end
%         end
        
%         function plot_grid(obj,variable_name,colocation)
%             if nargin<2
%                 variable_name=[];
%             end
%             if nargin<3
%                 colocation=0;
%             end
%             if isempty(variable_name)==0
%                 vindex=obj.stem_varset.get_index(variable_name);
%                 if isempty(vindex)==0
%                     figure;
%                     mapshow(obj.shape);
%                     hold on
%                     for i=1:obj.stem_varset.dim(vindex)
%                         plot(obj.stem_gridlist.grid{vindex}.coordinate(i,2),obj.stem_gridlist.grid{vindex}.coordinate(i,1),'*r');
%                     end
%                     box=obj.stem_gridlist.box;
%                     a=(box(4)-box(3))*0.2;
%                     b=(box(2)-box(1))*0.2;
%                     set(gca,'Xlim',[box(3)-a,box(4)+a]);
%                     set(gca,'Ylim',[box(1)-b,box(2)+b]);
%                     xlabel('Longitude');
%                     ylabel('Latitude');
%                     title([obj.stem_varset.name{vindex},' sites (',num2str(obj.stem_varset.dim(vindex)),')']);
%                 else
%                     error('The variable name is incorrect');
%                 end
%             else
%                 if colocation==0
%                     if obj.stem_varset.nvar<=3
%                         rows=1;
%                         cols=obj.stem_varset.nvar;
%                     else
%                         l=obj.stem_varset.nvar^0.5;
%                         if round(l)^2==obj.stem_varset.nvar
%                             rows=l;
%                             cols=l;
%                         else
%                             rows=ceil(l);
%                             cols=round(l);
%                         end
%                     end
%                     for i=1:obj.stem_varset.nvar
%                         subplot(rows,cols,i);
%                         %mapshow(obj.shape);%,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[.4 .4 .4]);
%                         %hold on
%                         %for j=1:length(obj.stem_gridlist.grid{i}.coordinate)
%                         %    plot(obj.stem_gridlist.grid{i}.coordinate(j,2),obj.stem_gridlist.grid{i}.coordinate(j,1),'*r');
%                         %end
%                         %box=obj.stem_gridlist.box;
%                         %a=(box(4)-box(3))*0.2;
%                         %b=(box(2)-box(1))*0.2;
%                         %set(gca,'Xlim',[box(3)-a,box(4)+a]);
%                         %set(gca,'Ylim',[box(1)-b,box(2)+b]);
% 
%                         axesm('MapProjection','sinusoid','MapLatLimit',[54.5 61],'MapLonLimit',[-8 0],...
%                             'MLineLocation',1,'PLineLocation',1,'GColor',[0 0 0],'ParallelLabel','on','MeridianLabel','on','Grid','on')
%                         obj.shape.X(isnan(obj.shape.X))=0;
%                         obj.shape.Y(isnan(obj.shape.Y))=0;
%                         indices=find(obj.shape.X==0);
%                         
%                         for j=1:length(indices)-1
%                             geoshow(obj.shape.Y(indices(j)+1:indices(j+1)-1),obj.shape.X(indices(j)+1:indices(j+1)-1),'DisplayType','polygon','FaceColor','none');
%                         end
%                         hold on
%                         geoshow(obj.stem_gridlist.grid{i}.coordinate(:,1),obj.stem_gridlist.grid{i}.coordinate(:,2),...
%                                 'DisplayType','multipoint','Marker','o','MarkerEdgeColor','black','MarkerFaceColor','black');
%                         title([obj.stem_varset.name{i},' (',num2str(obj.stem_varset.dim(i)),' sites)'],'FontSize',13);
%                         %xlabel('Longitude','FontSize',12);
%                         %ylabel('Latitude','FontSize',12);
%                     end
%                 else
%                     if obj.stem_varset.nvar>8
%                         error('Too many variables for this plot option');
%                     end
%                     colors='brygcmkw';
%                     scrsz = get(0,'ScreenSize');
%                     scrsz = scrsz + [100 100 -200 -200];
%                     figure('Position',scrsz);
%                     mapshow(obj.shape);%,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[.4 .4 .4]);
%                     hold on     
%                     for i=1:obj.stem_varset.nvar
%                         for j=1:length(obj.stem_gridlist.grid{i}.coordinate)
%                             h(i)=plot(obj.stem_gridlist.grid{i}.coordinate(j,2),obj.stem_gridlist.grid{i}.coordinate(j,1),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors(i),'MarkerSize',(obj.stem_varset.nvar-i+1)*4+3);
%                         end
%                     end
%                     box=obj.stem_gridlist.box;
%                     a=(box(4)-box(3))*0.2;
%                     b=(box(2)-box(1))*0.2;
%                     set(gca,'Xlim',[box(3)-a,box(4)+a]);
%                     set(gca,'Ylim',[box(1)-b,box(2)+b]);
%                     %title([obj.variable_name{i},' sites (',num2str(obj.var_dims(i)),')']);
%                     xlabel('Longitude');
%                     ylabel('Latitude');
%                     legend(h,obj.stem_varset.name,'Location','EastOutside');
%                 end
%             end
%         end
        
%         function hist_cov(obj,variable_name)
%             % display the histogram of the covariates for the specified variable
%             % variable_name: variable name string
%             vindex=obj.stem_varset.get_index(variable_name);
%             if not(isempty(vindex))
%                 l=obj.stem_covset.nbeta(vindex)^0.5;
%                 if round(l)^2==obj.stem_covset.nbeta(vindex)
%                     rows=l;
%                     cols=l;
%                 else
%                     rows=ceil(l);
%                     cols=round(l);
%                 end
%                 for i=1:obj.stem_covset.nbeta(vindex)
%                     subplot(rows,cols,i);
%                     a=obj.stem_covset.X{vindex}(:,i,:);
%                     hist(a(:));
%                     xlabel(obj.stem_covset.name{vindex}{(i)});
%                 end
%             else
%                 error('Variable name not recognized');
%             end            
%         end
        
%         function hist_var(obj)
%             % display the histogram of all the variables
%             if obj.stem_varset.nvar<=3
%                 rows=1;
%                 cols=obj.stem_varset.nvar;
%             else
%                 l=obj.stem_varset.nvar^0.5;
%                 if round(l)^2==obj.stem_varset.nvar
%                     rows=l;
%                     cols=l;
%                 else
%                     rows=ceil(l);
%                     cols=round(l);
%                 end
%             end
%             for i=1:obj.stem_varset.nvar
%                 subplot(rows,cols,i);
%                 hist(obj.stem_varset.Y{i}(:),20);
%                 xlabel(obj.stem_varset.name{i});
%             end
%         end
        
%         function temporal_cross_plot(obj,variable_name1,variable_name2,site,estimated)
%             if nargin<3
%                 error('Not enough input parameters');
%             end
%             var1_index=obj.stem_varset.get_index(variable_name1);
%             var2_index=obj.stem_varset.get_index(variable_name2);
%             if isempty(var1_index)
%                 error('The first variable_name is incorrect');
%             end
%             if isempty(var2_index)
%                 error('The second variable_name is incorrect');
%             end
%             if (nargin<4)||(isempty(site))
%                 site=1:length(obj.stem_gridlist.grid{var1_index}.coordinate);
%             end
%             if nargin<5
%                 estimated=0;
%             end
%             if estimated&&(isempty(obj.Y_hat))
%                 error('The model has not been estimated yet');
%             end
%             blocks=[0 cumsum(obj.stem_varset.dim)];
%             figure
%             done=0;
%             for j=site
%                 latlon1=obj.stem_gridlist.grid{var1_index}.coordinate(j,:);
%                 for i=1:length(obj.stem_gridlist.grid{var2_index}.coordinate)
%                     latlon2=obj.stem_gridlist.grid{var2_index}.coordinate(i,:);
%                     if sum(latlon1-latlon2)==0
%                         done=done+1;
%                         index1=blocks(var1_index)+j;
%                         index2=blocks(var2_index)+i;
%                         if estimated==0
%                             data1=obj.Y(index1,:);
%                             data2=obj.Y(index2,:);
%                         else
%                             data1=obj.Y_hat(index1,:);
%                             data2=obj.Y_hat(index2,:);
%                         end
%                         for t=1:length(data1)
%                             plot(data1(t),data2(t),'o','MarkerFaceColor',[t/length(data1) t/length(data1) t/length(data1)]);
%                             hold on
%                         end
%                     end
%                 end
%             end
%             if length(site)==1
%                 title(['Cross plot between ',variable_name1,' and ',variable_name2,' at site number',num2str(site)]);
%             else
%                 title(['Cross plot between ',variable_name1,' and ',variable_name2,' at ',num2str(done),' sites']);
%             end
%             xlabel(variable_name1);
%             ylabel(variable_name2);
%             if not(done)
%                 error('The two variables are not observed at the same site');
%             end
%         end
        
%         function v = temporal_cross_correlation(obj,type,graph)
%             % temporal cross-correlation between variables at sites. The variable with the lower number of sites is the reference variable
%             % graph: 0: no graph is plotted; 1: time-series graph and correlation graph are plotted
%             if nargin<2
%                 type='nearest';
%             end
%             if nargin<3
%                 graph=0;
%             end
%             if not(strcmp(type,'nearest')||strcmp(type,'colocated'))
%                 error('type must be either nearest or colocated');
%             end
%             v_matrix=eye(obj.stem_varset.nvar);
%             for j=1:obj.stem_varset.nvar-1
%                 for i=j+1:obj.stem_varset.nvar
%                     if obj.stem_varset.dim(i)<obj.stem_varset.dim(j)
%                         h=i;
%                         k=j;
%                         swtc=false;
%                     else
%                         h=j;
%                         k=i;
%                         swtc=true;
%                     end
%                     if graph
%                         figure
%                         l=obj.stem_varset.dim(h)^0.5;
%                         if round(l)^2==obj.stem_varset.dim(h)
%                             rows=l;
%                             cols=l;
%                         else
%                             rows=ceil(l);
%                             cols=round(l);
%                         end
%                     end
%                     ok_count=0;
%                     d=[];
%                     correlation=[];
%                     weight=[];
%                     for z=1:obj.stem_varset.dim(h)
%                         d(:,1)=abs(obj.stem_gridlist.grid{k}.coordinate(:,1)-obj.stem_gridlist.grid{h}.coordinate(z,1));
%                         d(:,2)=abs(obj.stem_gridlist.grid{k}.coordinate(:,2)-obj.stem_gridlist.grid{h}.coordinate(z,2));
%                         d=sum(d,2);
%                         [m,index]=min(d);
%                         if (strcmp(type,'colocated'))&&(sum(m)~=0)
%                             ok=false;
%                         else
%                             ok=true;
%                         end
%                         ok=true;
%                         if ok
%                             ok_count=ok_count+1;
%                             temp1=obj.stem_varset.Y{h}(z,:);
%                             temp2=obj.stem_varset.Y{k}(index,:);
%                             c=nancov(temp1,temp2);
%                             correlation(ok_count)=c(1,2)/(c(1,1)*c(2,2))^0.5;
%                             weight(ok_count)=sum(isnotnan(temp1)&isnotnan(temp2))/obj.stem_varset.T;
%                             if graph
%                                 subplot(rows,cols,z);
%                                 plot((temp1-nanmean(temp1))/nanstd(temp1),'b*');
%                                 hold on
%                                 plot((temp2-nanmean(temp2))/nanstd(temp2),'r*');
%                                 title(['Correlation: ',num2str(correlation(z))]);
%                             end
%                         end
%                     end
% %                     if graph
% %                         legend(obj.variable_name{h},obj.variable_name{k},'Location','EastOutside');
% %                     end
%                     L=isnotnan(correlation);
%                     correlation=correlation(L);
%                     weight=weight(L);
%                     if graph
%                         figure
%                         subplot(1,2,1);
%                         hist(correlation);
%                         xlabel('Correlation');
%                         subplot(1,2,2);
%                         plot((1-weight)*100,correlation,'*');
%                         xlabel('% full time-period missing');
%                         ylabel('Correlation');
%                     end
%                     mean_correlation=(correlation*weight')/sum(weight);
%                     v_matrix(h,k)=mean_correlation;
%                     v_matrix(k,h)=mean_correlation;
%                     if swtc
%                         disp(['Temporal mean weighted correlation between ',obj.stem_varset.name{h},' and ',obj.stem_varset.name{k},': ',num2str(mean_correlation), '  (',num2str(ok_count),' sites)']);
%                     else
%                         disp(['Temporal mean weighted correlation between ',obj.stem_varset.name{k},' and ',obj.stem_varset.name{h},': ',num2str(mean_correlation), '  (',num2str(ok_count),' sites)']);
%                     end
%                 end
%             end
%             if nargout>0
%                 v=v_matrix;
%                 if min(eig(v))<0
%                     warning('The output correlation matrix is not semi-definite positive!')
%                 end
%             end
%         end
        
%         function [matrix,average] = temporal_correlation(obj,variable_name)
%             var_index=obj.stem_varset.get_index(variable_name);    
%             if isempty(var_index)
%                 error('The variable name is incorrect')
%             end
%             data=obj.stem_varset.Y{var_index};
%             matrix=eye(size(data,1));
%             average=0;
%             counter=0;
%             for j=1:size(data,1)
%                 for i=j+1:size(data,1)
%                     a=data(i,:);
%                     b=data(j,:);
%                     L1=isnotnan(a);
%                     L2=isnotnan(b);
%                     L=L1&L2;
%                     a=a(L);
%                     b=b(L);
%                     temp=cov(a,b);
%                     matrix(i,j)=temp(1,2)/(var(a)*var(b))^0.5;
%                     matrix(j,i)=matrix(i,j);
%                     if isnotnan(matrix(i,j))
%                         average=average+matrix(i,j)*sum(L);
%                         counter=counter+sum(L);
%                     end
%                 end
%             end
%             average=average/counter;
%         end
        
%         function plot(obj,variable_name,time_step,plot_type,google_maps,covariate_name)
%             % display the animation for the selected variable or covariate
%             % variable_name: variable name string
%             % covariate_name: covariate name string. If provided the animation is on the covariate with respect to the specified
%             %                 variable; if not provided the animation is on the specified variable
%             % step_start: starting temporal step
%             % step_end: ending temportal step
%             if nargin<5
%                 error('Not enough input arguments');
%             end
%             if nargin<4
%                 plot_type='data';
%             end
%             if nargin<5
%                 google_maps=0;
%             end
%             if nargin<6
%                 covariate_name=[];
%             end
%             index_var=obj.stem_varset.get_index(variable_name);
%             if isempty(index_var)
%                 error('The variable name is incorrect');
%             end
%             if not(isempty(covariate_name))
%                 index_cov=obj.stem_covset.get_index(index_var,covariate_name);
%                 if isempty(index_cov)
%                     error('The covariate name is incorrect');
%                 end
%                 plot_cov=1;
%             else
%                 plot_cov=0;
%             end
%             if not(length(time_step)==2)
%                 error('time_step must be a 2x1 vector or cell array of date strings');
%             end
%             if iscell(time_step)
%                 date_num1=datenum(time_step{1},'dd-mm-yyyy');
%                 date_num2=datenum(time_step{2},'dd-mm-yyyy');
%                 tstep1=find(obj.stem_datestamp.stamp==date_num1);
%                 tstep2=find(obj.stem_datestamp.stamp==date_num2);
%                 if isempty(tstep1)
%                     error('Dates out of bounds');
%                 end
%                 if isempty(tstep2)
%                     error('Dates out of bounds');
%                 end
%                 if tstep2<tstep1
%                     error('Swap dates');
%                 end
%             else
%                 tstep1=time_step(1);
%                 tstep2=time_step(2);
%             end
%             if plot_cov
%                 error('Not yet supported for covariates');
%             else
%                 if strcmp(obj.stem_gridlist.grid{index_var}.grid_type,'regular')
%                     lat=reshape(obj.stem_gridlist.grid{index_var}.coordinate(:,1),obj.stem_gridlist.grid{index_var}.grid_size);
%                     long=reshape(obj.stem_gridlist.grid{index_var}.coordinate(:,2),obj.stem_gridlist.grid{index_var}.grid_size);
%                     data=obj.stem_varset.Y{index_var}(:,tstep);
%                     data=reshape(data,obj.stem_gridlist.grid{index_var}.grid_size);
%                     mapshow(obj.shape);
%                     hold on
%                     mapshow(long,lat,data,'DisplayType','texture');
%                     set(gca,'Xlim',[obj.stem_gridlist.box(3),obj.stem_gridlist.box(4)]);
%                     set(gca,'Ylim',[obj.stem_gridlist.box(1),obj.stem_gridlist.box(2)]);
%                     title(datestr(obj.stem_datestamp.stamp(t)));
%                 else
%                     full_data=obj.stem_varset.Y{index_var};
%                     if strcmp(plot_type,'data')
%                         data=full_data;
%                     end
%                     if strcmp(plot_type,'average');
%                         data=nanmean(full_data,2);
%                         tstep1=1;
%                         tstep2=1;
%                     end
%                     if strcmp(plot_type,'std');
%                         data=nanstd(full_data,0,2);
%                         tstep1=1;
%                         tstep2=1;
%                     end
%                     if strcmp(plot_type,'autocorr')
%                         for i=1:size(full_data,1)
%                             temp=full_data(i,:);
%                             temp=temp(isnotnan(temp));
%                             xc=xcorr(temp,'coeff');
%                             l=length(xc);
%                             xc=xc((l+1)/2:end);
%                             data(i,1)=sum(xc(2:2+5));
%                         end
%                         tstep1=1;
%                         tstep2=1;
%                     end
%                     min_value=min(data(:));
%                     max_value=max(data(:));
%                     if not(google_maps)
%                         for t=tstep1:tstep2
%                             hold off
%                             mapshow(obj.shape);
%                             hold on
%                             for i=1:length(obj.stem_gridlist.grid{index_var}.coordinate)
%                                 if isnotnan(data(i,t))
%                                     value=(data(i,t)-min_value)/(max_value-min_value)*15+5;
%                                     plot(obj.stem_gridlist.grid{index_var}.coordinate(i,2),obj.stem_gridlist.grid{index_var}.coordinate(i,1),'ro','MarkerSize',round(value));
%                                 end
%                             end
%                             set(gca,'Xlim',[obj.stem_gridlist.box(3),obj.stem_gridlist.box(4)]);
%                             set(gca,'Ylim',[obj.stem_gridlist.box(1),obj.stem_gridlist.box(2)]);
%                             grid on
%                             if strcmp(plot_type,'data')
%                                 title([variable_name,' for ',datestr(obj.stem_datestamp.stamp(t))]);
%                             end
%                             if strcmp(plot_type,'average')
%                                 title(['Time average for ',variable_name]);
%                             end
%                             if strcmp(plot_type,'std')
%                                 title(['Standard deviation for ',variable_name]);
%                             end
%                         end
%                     else
%                         datafile=[tstep2-tstep1+1,obj.stem_datestamp.stamp(tstep1)];
%                         csvwrite('..\Data\google_bridge\parameters.csv',datafile);
%                         for t=tstep1:tstep2
%                             value=(data(:,t)-min_value)/(max_value-min_value)*5+1;
%                             datafile=[obj.stem_gridlist.grid{index_var}.coordinate(:,1),obj.stem_gridlist.grid{index_var}.coordinate(:,2),value,data(:,t)];
%                             csvwrite(['..\Data\google_bridge\stations',num2str(t-tstep1+1),'.csv'],datafile);
%                         end
%                         winopen('..\Data\google_bridge\open.bat');
%                     end
%                 end
%             end           
%         end
        
%       if color_value<0.5
%         color=[color_value*2 1 0];
%       else
%         color=[1 1-(color_value-0.5)*2 0];
%       end
    end
end

