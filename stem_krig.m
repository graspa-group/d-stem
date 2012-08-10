%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_krig < handle
    %stem kriging class
    
    properties
        stem_model=[];
        X_all=[];
        idx_notnan=[];
    end
    
    methods
        
        function obj = stem_krig(stem_model)
            if nargin<1
                error('All the input arguments must be provided');
            end
            obj.stem_model=stem_model;
        end
        
        function st_krig_result = kriging(obj,block_size,variable_name,grid,mask,X,back_transform,no_varcov,crossval)
            %This is the function that actually implements the kriging
            
            %block_size: when the number of kriging sites is large the
            %kriging can be implemented block by block with dimension given
            %by block_size>=0. If block_size=0 all the kriging sites are
            %considered in one block.
            
            %variable_name: the name of the variable to be kriged
            
            %grid: the stem_grid object with the sites where the kriging
            %has to be implemented
            
            %mask: a column vector with dimension the number of sites in
            %grid the elements of which are NaN if kriging is not needed at
            %the respective site
            
            %X: the kriging covariates. X has to be provided when the model
            %has space-time varying coefficients (covariates). X can be
            %either a structure or the name of a folder. When X is a
            %structure it must contain:
                %X_krig.X_all: the matrix of all the coefficients needed
                %X_krig.name: the names of each column of coefficients
                %X_krig.date_stamp: a stem_datestamp object
            %When X is a folder name it must contain the name of the folder
            %where the coefficients are saved block by block. The blocks
            %must contain only the coefficients related to the non-masked
            %sites. The structure of the blocks is specified in ????
            
            %back_transform: is a flag that specify if the kriged variable
            %has to be back transformed
            
            %no_varcov: when this flag is 1, the kriging
            %variance-covariance is not evaluated (to save time)
            
            %crossval: when this flag is 1, the variable is kriged over the
            %cross-validation sites. In this case X and grid should not be
            %provided. NOTE THAT THE CROSS-VALIDATION KRIGING IS
            %AUTOMATICALLY IMPLEMENTED AFTER THE EM-ESTIMATION IN THE CASE
            %THE CROSS-VALIDATION INFORMATION IS PROVIDED!
            
            if nargin<4
                error('Not enough input arguments');
            end
            if isempty(obj.stem_model)
                error('The stem_model property is not setted');
            end
            if block_size<0
                error('block_size cannot be < 0');
            end
            index_var=obj.stem_model.stem_data.stem_varset_g.get_Y_index(variable_name);
            if isempty(index_var)
                error('The variable name is incorrect');
            end
            if nargin<9
                crossval=0;
            end
            if (crossval)&&not(obj.stem_model.cross_validation)
                error('The stem_model object does not contain cross-validation information');
            end
            if (crossval)&&(not(isempty(X)))
                disp('WARNING: the X provided in not considered as the covariates of cross validation are used');
            end
            if (crossval)&&(not(isempty(grid)))
                disp('WARNING: the grid provided in not considered as the grid of cross validation is used');
            end                
            if crossval
                grid=obj.stem_model.stem_data.stem_crossval.stem_gridlist.grid{1};
            end
            
            if nargin<5
                mask=[];
                obj.idx_notnan=[1:length(grid.coordinate)]';
            else
                if not(isempty(mask))
                    obj.idx_notnan=find(isnotnan(mask(:)));
                else
                    obj.idx_notnan=[1:length(grid.coordinate)]';
                end
            end            
            if nargin<6
                X_all=[];
            end
            if (nargin<7)||isempty(back_transform)
                back_transform=1;
            end
       
            if block_size==0
                blocks_krig=[0 size(obj.idx_notnan,1)];
            else
                blocks_krig=0:block_size:size(obj.idx_notnan,1);
                if not(blocks_krig(end)==size(obj.idx_notnan,1))
                    blocks_krig=[blocks_krig size(obj.idx_notnan,1)];
                end
            end
            
            if not(crossval)
                %indexes recovering and block test
                loadfromfile=0;
                if not(isempty(X))
                    if isstruct(X)
                        %the X are directly provided
                        if not(isempty(mask))
                            if not(size(mask,1)*size(mask,2)==size(grid.coordinate,1))
                                error('The mask size does not match the number of coordinates');
                            end
                        end
                        if not(isempty(X.X_all))
                            if not(size(X.X_all,1)==size(grid.coordinate,1))
                                error('The size of X_all does not match the number of coordinates');
                            end
                            obj.X_all=X.X_all(obj.idx_notnan,:,:);
                        end
                        
                        idx_datestamp=[];
                        for j=1:length(obj.stem_model.stem_data.stem_datestamp.stamp)
                            idx=find(X.date_stamp.stamp==obj.stem_model.stem_data.stem_datestamp.stamp(j),1);
                            if isempty(idx)
                                error('The kriging block does not contain the correct datestamp');
                            end
                            idx_datestamp=[idx_datestamp,idx];
                        end
                        
                        %recover the indices
                        if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_rg_name))
                            X_rg_name=obj.stem_model.stem_data.stem_varset_g.X_rg_name{index_var};
                            if not(isempty(X_rg_name))
                                idx_rg=[];
                                for j=1:length(X_rg_name)
                                    cmp=strcmp(X.name,X_rg_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_rg=[idx_rg,find(cmp,1)];
                                    end
                                end
                            else
                                idx_rg=[];
                            end
                        else
                            idx_rg=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_g_name))
                            X_g_name=obj.stem_model.stem_data.stem_varset_g.X_g_name{index_var};
                            if not(isempty(X_g_name))
                                idx_g=[];
                                for j=1:length(X_g_name)
                                    cmp=strcmp(X.name,X_g_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_g=[idx_g,find(cmp,1)];
                                    end
                                end
                            else
                                idx_g=[];
                            end
                        else
                            idx_g=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_beta_name))
                            X_beta_name=obj.stem_model.stem_data.stem_varset_g.X_beta_name{index_var};
                            if not(isempty(X_beta_name))
                                idx_beta=[];
                                for j=1:length(X_beta_name)
                                    cmp=strcmp(X.name,X_beta_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_beta=[idx_beta,find(cmp,1)];
                                    end
                                end
                            else
                                idx_beta=[];
                            end
                        else
                            idx_beta=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_time_name))
                            X_time_name=obj.stem_model.stem_data.stem_varset_g.X_time_name{index_var};
                            if not(isempty(X_time_name))
                                idx_time=[];
                                for j=1:length(X_time_name)
                                    cmp=strcmp(X.name,X_time_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_time=[idx_time,find(cmp,1)];
                                    end
                                end
                            else
                                idx_time=[];
                            end
                        else
                            idx_time=[];
                        end
                    else
                        %a folder name has been provided
                        loadfromfile=1;
                        folder=X;
                        if block_size==0
                            error('A block size >0 must be provided');
                        end
                        files=dir([X,'*.mat']);
                        %test data compatibility
                        disp('Hard disk kriging blocks evaluation started...');
                        counter=0;
                        for i=1:length(files)
                            load([X,files(i).name]); %load X_krig_block variable
                            idx_datestamp=[];
                            for j=1:length(obj.stem_model.stem_data.stem_datestamp.stamp)
                                idx=find(block.date_stamp==obj.stem_model.stem_data.stem_datestamp.stamp(j),1);
                                if isempty(idx)
                                    error('The kriging block does not contain the correct datestamp');
                                end
                                idx_datestamp=[idx_datestamp,idx];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_rg_name))
                                X_rg_name=obj.stem_model.stem_data.stem_varset_g.X_rg_name{index_var};
                                if not(isempty(X_rg_name))
                                    idx_rg=[];
                                    for j=1:length(X_rg_name)
                                        cmp=strcmp(block.label,X_rg_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_rg=[idx_rg,find(cmp,1)];
                                        end
                                    end
                                else
                                    idx_rg=[];
                                end
                            else
                                idx_rg=[];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_g_name))
                                X_g_name=obj.stem_model.stem_data.stem_varset_g.X_g_name{index_var};
                                if not(isempty(X_g_name))
                                    idx_g=[];
                                    for j=1:length(X_g_name)
                                        cmp=strcmp(block.label,X_g_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_g=[idx_g,find(cmp,1)];
                                        end
                                    end
                                else
                                    idx_g=[];
                                end
                            else
                                idx_g=[];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_beta_name))
                                X_beta_name=obj.stem_model.stem_data.stem_varset_g.X_beta_name{index_var};
                                if not(isempty(X_beta_name))
                                    idx_beta=[];
                                    for j=1:length(X_beta_name)
                                        cmp=strcmp(block.label,X_beta_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_beta=[idx_beta,find(cmp,1)];
                                        end
                                    end
                                else
                                    idx_beta=[];
                                end
                            else
                                idx_beta=[];
                            end
                            
                            if not(isempty())
                                X_time_name=obj.stem_model.stem_data.stem_varset_g.X_time_name{index_var};
                                if not(isempty(X_time_name))
                                    idx_time=[];
                                    for j=1:length(X_time_name)
                                        cmp=strcmp(block.label,X_time_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_time=[idx_time,find(cmp,1)];
                                        end
                                    end
                                else
                                    idx_time=[];
                                end
                            else
                                idx_time=[];
                            end
                            
                            if not(size(block.data,1)==blocks_krig(i+1)-blocks_krig(i))
                                error(['Wrong number of rows in kriging block',num2str(i)]);
                            end
                            counter=counter+size(block.data,1);
                        end
                        if not(counter==size(obj.idx_notnan,1))
                            error('The total number of rows in kriging blocks is wrong');
                        end
                        disp('Hard disk kriging blocks evaluation ended');
                    end
                end
            else
                %kriging using cross-validation data
            end
            
            
            disp('Kriging started...');
            K=obj.stem_model.stem_par.k;
            
            st_krig_result=stem_krig_result(variable_name,grid,obj.stem_model.stem_data.shape);
            st_krig_result.y_hat=zeros(size(grid.coordinate,1),obj.stem_model.T);
            if not(no_varcov)
                st_krig_result.var_y_hat=zeros(size(grid.coordinate,1),obj.stem_model.T);
            end
            if K>0
                st_krig_result.E_wg_y1=zeros(size(grid.coordinate,1),obj.stem_model.T,K);
            end
            
            for i=1:length(blocks_krig)-1
                ct1=clock;
                disp(['Kriging block ',num2str(i),' of ',num2str(length(blocks_krig)-1)]);
                block_krig=(blocks_krig(i)+1):blocks_krig(i+1);
                block_krig_length=length(block_krig);
                blocks=cumsum(obj.stem_model.stem_data.stem_varset_g.dim);
                
                %Y manage
                Y_add=nan(block_krig_length,obj.stem_model.T);
                obj.stem_model.stem_data.stem_varset_g.Y{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.Y{index_var},Y_add);
                clear Y_add
                
                %X manage
                if not(crossval)
                    if not(isempty(obj.X_all))||loadfromfile
                        if loadfromfile
                            load([folder,files(i).name]);
                            block.data=block.data(:,:,idx_datestamp);
                        else
                            if size(obj.X_all,3)>1
                                block.data=obj.X_all(block_krig,:,idx_datestamp);
                            else
                                block.data=obj.X_all(block_krig,:,1);
                            end
                            block.lat=grid.coordinate(obj.idx_notnan(block_krig),1);
                            block.lon=grid.coordinate(obj.idx_notnan(block_krig),2);
                        end
                        
                        if not(isempty(idx_rg))
                            X_krig_block=block.data(:,idx_rg,:);
                            if obj.stem_model.stem_data.stem_varset_g.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_g.X_rg_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_g.X_rg_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_g.X_rg{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                                warning('Only the first temporal step is considered');
                            end
                            obj.stem_model.stem_data.stem_varset_g.X_rg{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_rg{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idx_g))
                            X_krig_block=block.data(:,idx_g,:);
                            if obj.stem_model.stem_data.stem_varset_g.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_g.X_g_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_g.X_g_stds{index_var}(j);
                                end
                            end
                            temp=[];
                            for k=1:size(X_krig_block,2)
                                temp1=X_krig_block(:,k,:);
                                temp=cat(4,temp,temp1);
                            end
                            X_krig_block=temp;
                            if size(obj.stem_model.stem_data.stem_varset_g.X_g{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1,:);
                                warning('Only the first temporal step is considered');
                            end
                            obj.stem_model.stem_data.stem_varset_g.X_g{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_g{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idx_beta))
                            X_krig_block=block.data(:,idx_beta,:);
                            if obj.stem_model.stem_data.stem_varset_g.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_g.X_beta_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_g.X_beta_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_g.X_beta{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                                warning('Only the first temporal step is considered');
                            end
                            obj.stem_model.stem_data.stem_varset_g.X_beta{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_beta{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idx_time))
                            X_krig_block=block.data(:,idx_time,:);
                            if obj.stem_model.stem_data.stem_varset_g.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_g.X_time_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_g.X_time_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_g.X_time{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                                warning('Only the first temporal step is considered');
                            end
                            obj.stem_model.stem_data.stem_varset_g.X_time{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_time{index_var},X_krig_block);
                        end
                    end
                else
                    %cross-validation data
                    block.lat=grid.coordinate(obj.idx_notnan(block_krig),1);
                    block.lon=grid.coordinate(obj.idx_notnan(block_krig),2);
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset.X_rg))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset.X_rg{index_var}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_g.X_rg{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_rg{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset.X_g))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset.X_g{index_var}(obj.idx_notnan(block_krig),:,:,:);
                        obj.stem_model.stem_data.stem_varset_g.X_g{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_g{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset.X_beta))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset.X_beta{index_var}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_g.X_beta{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_beta{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset.X_time))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset.X_time{index_var}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_g.X_time{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_g.X_time{index_var},X_krig_block);
                    end
                end
                
                %Grid manage
                obj.stem_model.stem_data.stem_gridlist_g.grid{index_var}.coordinate=cat(1,obj.stem_model.stem_data.stem_gridlist_g.grid{index_var}.coordinate,[block.lat,block.lon]);
                obj.stem_model.stem_data.update_data();
                obj.stem_model.stem_data.update_distance('ground');
                if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                    obj.stem_model.stem_data.update_M();
                end
                
                %kriging
                [y_hat,var_y_hat,E_wg_y1]=obj.E_step();
                st_krig_result.y_hat(obj.idx_notnan(block_krig),:)=y_hat(blocks(index_var)+1:blocks(index_var)+block_krig_length,:);
                if not(no_varcov)
                    st_krig_result.var_y_hat(obj.idx_notnan(block_krig),:)=var_y_hat(blocks(index_var)+1:blocks(index_var)+block_krig_length,:);
                end
                if K>0
                    st_krig_result.E_wg_y1(obj.idx_notnan(block_krig),:,:)=E_wg_y1(blocks(index_var)+1:blocks(index_var)+block_krig_length,:,:);
                end
                
                %restore original
                obj.stem_model.stem_data.stem_varset_g.Y{index_var}(end-block_krig_length+1:end,:)=[];
                if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_rg))
                    obj.stem_model.stem_data.stem_varset_g.X_rg{index_var}(end-block_krig_length+1:end,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_g))
                    obj.stem_model.stem_data.stem_varset_g.X_g{index_var}(end-block_krig_length+1:end,:,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_beta))
                    obj.stem_model.stem_data.stem_varset_g.X_beta{index_var}(end-block_krig_length+1:end,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_g.X_time))
                    obj.stem_model.stem_data.stem_varset_g.X_time{index_var}(end-block_krig_length+1:end,:,:)=[];
                end
                obj.stem_model.stem_data.stem_gridlist_g.grid{index_var}.coordinate(end-block_krig_length+1:end,:)=[];
                obj.stem_model.stem_data.update_data(); %verificare se conviene salvare invece di ricalcolare
                obj.stem_model.stem_data.update_distance(); %verificare se conviene salvare invece di ricalcolare
                ct2=clock;
                disp(['Kriging block ended in ',stem_time(etime(ct2,ct1))]);
            end


            if back_transform&&obj.stem_model.stem_data.stem_varset_g.standardized
                disp('Back-transformation...');
                s=obj.stem_model.stem_data.stem_varset_g.Y_stds{index_var};
                m=obj.stem_model.stem_data.stem_varset_g.Y_means{index_var};
                if (obj.stem_model.stem_data.stem_varset_g.standardized)&&not(obj.stem_model.stem_data.stem_varset_g.log_transformed)
                    st_krig_result.y_hat=st_krig_result.y_hat*s+m;
                    if not(no_varcov)
                        st_krig_result.var_y_hat=st_krig_result.var_y_hat*s^2;
                    end
                end
                if (obj.stem_model.stem_data.stem_varset_g.standardized)&&(obj.stem_model.stem_data.stem_varset_g.log_transformed)
                    y_hat=st_krig_result.y_hat;
                    var_y_hat=st_krig_result.var_y_hat;
                    st_krig_result.y_hat=exp(y_hat*s+m+(s^2)/2);
                    if not(no_varcov)
                        st_krig_result.var_y_hat=(var_y_hat*s^2)*(exp(m)^2);
                    end
                end
                disp('Back-transformation ended.');
            end
            
            if strcmp(grid.grid_type,'regular')
                disp('Data reshaping...');
                st_krig_result.y_hat=reshape(st_krig_result.y_hat,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T);
                if not(no_varcov)
                    st_krig_result.var_y_hat=reshape(st_krig_result.var_y_hat,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T);
                end
                if K>0
                    st_krig_result.E_wg_y1=reshape(st_krig_result.E_wg_y1,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T,K);
                end
                disp('Date reshaping ended.');
            end
            
            
            if not(isempty(mask))&&strcmp(grid.grid_type,'regular')
                disp('Applying mask...');
                mask(isnotnan(mask))=1;
                for t=1:size(st_krig_result.y_hat,3)
                    st_krig_result.y_hat(:,:,t)=st_krig_result.y_hat(:,:,t).*mask;
                    if not(no_varcov)
                        st_krig_result.var_y_hat(:,:,t)=st_krig_result.var_y_hat(:,:,t).*mask;
                    end
                    for k=1:K
                        st_krig_result.E_wg_y1(:,:,t,k)=st_krig_result.E_wg_y1(:,:,t,k).*mask;
                    end
                end
                disp('Mask applied.');
            end
            
            st_krig_result.variable_name=variable_name;
        end
        
        function [y_hat,var_y_hat,E_wg_y1] = E_step(obj)
            
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                Nr=obj.stem_model.stem_data.stem_varset_r.N;
            else
                Nr=0;
            end
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            T=obj.stem_model.stem_data.T;
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            if obj.stem_model.tapering
                %in the case of tapering system_size is increased to avoid the computation of only the diagonal
                obj.stem_model.set_system_size(1e10);
            end
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            if p>0
                st_kalmansmoother_result=obj.stem_model.stem_EM_result.stem_kalmansmoother_result;
                if not(data.X_time_tv)
                    if obj.stem_model.tapering
                        var_Zt=sparse(data.X_time(:,:,1))*sigma_Z*sparse(data.X_time(:,:,1)'); 
                    else
                        var_Zt=data.X_time(:,:,1)*sigma_Z*data.X_time(:,:,1)'; 
                    end
                end
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo+var_Zt;
                end                
            else
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[]);    
                var_Zt=[];
                %variance of Y
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end                
            end
           
            res=data.Y;
            if not(isempty(data.X_beta))
                Xbeta=zeros(N,T);
                if data.X_beta_tv
                    for t=1:T
                        Xbeta(:,t)=data.X_beta(:,:,t)*par.beta;
                    end
                else
                    Xbeta=repmat(data.X_beta(:,:,1)*par.beta,1,T);
                end
                res=res-Xbeta;
                y_hat=Xbeta;
            else
                y_hat=zeros(size(res));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            %cov_wr_yz time invariant case
            if not(isempty(data.X_rg))
                if not(data.X_rg_tv)
                    cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,1),'r'),aj_rg,'r');
                end
                E_wr_y1=zeros(Nr,T);
                diag_Var_wr_y1=zeros(Nr,T);
                cov_wr_z_y1=zeros(Nr,p,T);
            end
            
            %cov_wg_yz time invariant case
            if not(isempty(data.X_g))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_g{1});
                    [Ig,Jg]=ind2sub(size(sigma_W_g{1}),Lg);
                end
                if not(data.X_g_tv)
                    for k=1:K
                        cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,1,k),'r'),aj_g(:,k),'r');
                    end
                end
                for h=1:K
                    for k=h+1:K
                        cov_wgk_wgh_y1{k,h}=zeros(Ng,T);
                    end
                end
                E_wg_y1=zeros(Ng,T,K);
                diag_Var_wg_y1=zeros(Ng,T,K);
                cov_wg_z_y1=zeros(Ng,p,T,K);
            else
                E_wg_y1=[];
            end
            
            if not(isempty(data.X_rg)) && not(isempty(data.X_g))
                M_cov_wr_wg_y1=zeros(N,T,K);
            else
                M_cov_wr_wg_y1=[];
            end
            
            diag_Var_e_y1=zeros(N,T);
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(data.Y(:,t)));
                
                if data.X_rg_tv
                    tRG=t;
                else
                    tRG=1;
                end
                if data.X_time_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_g_tv
                    tG=t;
                else
                    tG=1;
                end                  
                
                %evaluate var_yt in the time variant case
                if data.X_tv
                    if not(isempty(data.X_rg))
                        sigma_geo=D_apply(D_apply(M_apply(sigma_W_r,M,'b'),data.X_rg(:,1,tRG),'b'),aj_rg,'b');
                    end
                    
                    if not(isempty(data.X_g))
                        if isempty(data.X_rg)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(data.X_g,4)
                            sigma_geo=sigma_geo+D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'b'),aj_g(:,k),'b');
                        end
                    end
                    if isempty(data.X_g)&&isempty(data.X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if not(isempty(data.X_time))
                        if data.X_time_tv
                            if obj.stem_model.tapering
                                var_Zt=sparse(data.X_time(:,:,tT))*sigma_Z*sparse(data.X_time(:,:,tT)');
                            else
                                var_Zt=data.X_time(:,:,tT)*sigma_Z*data.X_time(:,:,tT)';
                            end
                        end
                        var_Yt=sigma_geo+var_Zt;
                    else
                        var_Yt=sigma_geo;
                    end
                end
                
                %check if the temporal loadings are time variant
                if not(isempty(data.X_time))
                    if obj.stem_model.tapering
                        temp=sparse(data.X_time(:,:,tT))*st_kalmansmoother_result.Pk_s(:,:,t+1);
                    else
                        temp=data.X_time(:,:,tT)*st_kalmansmoother_result.Pk_s(:,:,t+1);
                    end
                    if N>obj.stem_model.system_size
                        for i=1:size(diag_Var_e_y1,1)
                            %update diag(Var(e|y1))
                            diag_Var_e_y1(i,t)=temp(i,:)*data.X_time(i,:,tT)';
                        end
                    else
                        if obj.stem_model.tapering
                            diag_Var_e_y1(:,t)=diag(sparse(temp)*sparse(data.X_time(:,:,tT)'));
                        else
                            diag_Var_e_y1(:,t)=diag(temp*data.X_time(:,:,tT)');
                        end
                    end
                end
                
                %build the Ht matrix
                if not(isempty(var_Zt))
                    temp=st_kalmansmoother_result.zk_s(:,t+1);
                    H1t=[var_Yt(Lt,Lt), data.X_time(Lt,:,tT)*sigma_Z; sigma_Z*data.X_time(Lt,:,tT)', sigma_Z];
                    y_hat(:,t)=y_hat(:,t)+data.X_time(:,:,tT)*temp;
                else
                    H1t=var_Yt(Lt,Lt);
                    temp=[];
                end
                
                if obj.stem_model.tapering
                    cs=[];
                    r = symamd(H1t);
                    chol_H1t=chol(H1t(r,r));
                    temp2=[res(Lt,t);temp];
                    cs(r,1)=chol_solve(chol_H1t,temp2(r));
                else
                    chol_H1t=chol(H1t);
                    cs=chol_solve(chol_H1t,[res(Lt,t);temp]);
                end
                
                if not(isempty(data.X_rg))
                    %check if the remote loadings are time variant
                    if data.X_rg_tv
                        %cov_wr_yz time variant case
                        cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,tRG),'r'),aj_rg,'r');
                    end
                    cov_wr_y1z=[cov_wr_y(:,Lt),zeros(size(cov_wr_y,1),p)];
                                      
                    %compute E(w_r|y1);
                    E_wr_y1(:,t)=cov_wr_y1z*cs;
                    %compute diag(Var(w_r|y1))
                    if obj.stem_model.tapering
                        temp_r(r,:)=chol_solve(chol_H1t,cov_wr_y1z(:,r)');
                    else
                        temp_r=chol_solve(chol_H1t,cov_wr_y1z');
                    end                    
                    
                    %DIVERSO RISPETTO A E-STEP DI STEM_EM IN QUANTO SERVE SOLO LA DIAGONALE!!
                    for i=1:size(temp_r,2)
                        diag_Var_wr_y1(i,t)=sigma_W_r(i,i)-cov_wr_y1z(i,:)*temp_r(:,i);
                    end
                    
                    if p>0
                        %compute cov(w_r,z|y1)
                        cov_wr_z_y1(:,:,t)=temp_r(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                        for i=1:size(temp_r,2)
                            diag_Var_wr_y1(i,t)=diag_Var_wr_y1(i,t)+cov_wr_z_y1(i,:,t)*temp_r(end-p+1:end,i);
                        end
                        %update diag(Var(e|y1))
                        temp=D_apply(D_apply(M_apply(cov_wr_z_y1(:,:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                        if N>obj.stem_model.system_size
                            for i=1:size(temp,1)
                                diag_Var_e_y1(i,t)=diag_Var_e_y1(i,t)+2*temp(i,:)*data.X_time(i,:,tT)'; %notare 2*
                            end
                        else
                            %faster for N small
                            if obj.stem_model.tapering
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(sparse(temp)*sparse(data.X_time(:,:,tT)'));
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)');
                            end
                        end
                    else
                        cov_wr_z_y1=[];
                    end
                    %update y_hat
                    y_hat(:,t)=y_hat(:,t)+D_apply(D_apply(M_apply(E_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(M_apply(diag_Var_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'b'),aj_rg,'b'); %tested
                end
                
                if not(isempty(data.X_g))
                    %check if the ground loadings are time variant
                    if data.X_g_tv    
                        %cov_wg_yz time invariant case
                        for k=1:K
                            cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'r'),aj_g(:,k),'r');
                        end
                    end
                    for k=1:K
                        cov_wg_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                        %compute E(w_g_k|y1);
                        E_wg_y1(:,t,k)=cov_wg_y1z*cs;
                        %compute diag(Var(w_g_k|y1))
                        if obj.stem_model.tapering
                            temp_g{k}(r,:)=chol_solve(chol_H1t,cov_wg_y1z(:,r)');
                        else
                            temp_g{k}=chol_solve(chol_H1t,cov_wg_y1z');
                        end
                        %DIVERSO RISPETTO A E-STEP DI STEM_EM IN QUANTO SERVE SOLO LA DIAGONALE!!
                        for i=1:size(temp_g{k},2)
                            %VERIFICARE SE E' NECESSARIO TRIMMARE TEMP_G{K} IN QUANTO DIVERSO DA E_STEP DI STEM_EM!!
                            diag_Var_wg_y1(:,t,k)=sigma_W_g{k}(i,i)-cov_wg_y1z(i,:)*temp_g{k}(:,i);
                        end

                        if p>0
                            %compute cov(w_g,z|y1)
                            cov_wg_z_y1(:,:,t,k)=temp_g{k}(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                            for i=1:size(temp_g,2)
                                diag_Var_wg_y1(:,t,k)=diag_Var_wg_y1(:,t,k)+cov_wg_z_y1(i,:,t,k)*temp_g{k}(end-p+1:end,i);
                            end
                            %update diag(Var(e|y1))
                            temp=D_apply(D_apply(cov_wg_z_y1(:,:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                            if N>obj.stem_model.system_size
                                for i=1:size(temp,1)
                                    diag_Var_e_y1(i,t)=diag_Var_e_y1(i,t)+2*temp(i,:)*data.X_time(i,:,tT)'; %notare 2*
                                end
                            else
                                if obj.stem_model.tapering
                                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(sparse(temp)*sparse(data.X_time(:,:,tT)')); 
                                else
                                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)'); 
                                end
                            end
                        else
                            cov_wg_z_y1=[];
                        end
                        %y_hat
                        y_hat(:,t)=y_hat(:,t)+D_apply(D_apply(E_wg_y1(:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(diag_Var_wg_y1(:,t,k),data.X_g(:,:,tG,k),'b'),aj_g(:,k),'b'); %K varianze                        
                        
                        if not(isempty(data.X_rg))
                            %compute M_cov(w_r,w_g|y1); cioè M*cov(w_r,w_g|y1) da tenere in considerazione nelle forme chiuse!
                            if length(M)>obj.stem_model.system_size
                                for i=1:length(M)
                                    %tested
                                    if p>0
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i)+cov_wr_z_y1(M(i),:,t)*temp_g{k}(end-p+1:end,i); %ha già l'M_apply su left!!
                                    else
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i);
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M))+cov_wr_z_y1(M,:,t)*temp_g{k}(end-p+1:end,1:length(M))); %ha già l'M_apply su left!!
                                else
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1)) - tested
                            temp=D_apply(D_apply(M_cov_wr_wg_y1(:,t,k),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                            temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,k);zeros(Nr,1)],'l'),aj_g(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_gk,w_gh|y1);
                        for h=1:K
                            for k=h+1:K
                                cov_wgk_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                                if N>obj.stem_model.system_size
                                    for i=1:size(cov_wgk_y1z,1)
                                        if not(isempty(cov_wg_z_y1))
                                            cov_wgk_wgh_y1{k,h}(i,t)=-cov_wgk_y1z(i,:)*temp_g{h}(:,i)+cov_wg_z_y1(i,:,t,k)*temp_g{h}(end-p+1:end,i);
                                        else
                                            cov_wgk_wgh_y1{k,h}(i,t)=-cov_wgk_y1z(i,:)*temp_g{h}(:,i);
                                        end
                                    end
                                else
                                    if not(isempty(cov_wg_z_y1))
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h}+cov_wg_z_y1(:,:,t,k)*temp_g{h}(end-p+1:end,:));
                                    else
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h});
                                    end
                                end
                                temp=D_apply(D_apply(cov_wgk_wgh_y1{k,h}(:,t),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                                temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,h);zeros(Nr,1)],'l'),aj_g(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    end
                end
                clear temp_g
            end
            var_y_hat=diag_Var_e_y1;
        end
           
%         function st_bootkrig_result = rndpar_kriging(obj,block_size,variable_name,X,sub,mask,back_transform,replications,directory,no_varcov)
%             if nargin<3
%                 error('The variable name to be kriged must be provided');
%             end
%             if nargin<4
%                 X=[];
%             end
%             if nargin<5
%                 sub=1;
%             end
%             if isempty(sub)==0
%                 sub=1;
%             end
%             if nargin<6
%                 mask=[];
%                 obj.idx_notnan=ones(length(obj.Grid),1);
%             else
%                 obj.idx_notnan=find(isnotnan(mask(:)));
%             end
%             if (nargin<7)||isempty(back_transform)
%                 back_transform=1;
%             end
%             if nargin<8
%                 replications=100;
%             end
%             
%             if (nargin>=9)&&not(isempty(directory))
%                 savetofile=1;
%                 if not(strcmp(directory(end),'/'))
%                     directory=[directory,'/'];
%                 end
%             else
%                 savetofile=0;
%             end
%             
%             if nargin<10
%                 no_varcov=0;
%             end
% 
%             if isempty(obj.stem_model)
%                 error('The stem_model property is not setted');
%             end
%             if not(isempty(obj.stem_model.stem_data.X))&&isempty(X)
%                 error('Covariates must be provided for this model');
%             end
%             
%             index=obj.stem_model.stem_data.stem_varset.get_index(variable_name);
%             if isempty(index)
%                 error('The variable name is incorrect');
%             end
%             
%             if block_size==0
%                 blocks_krig=[0 size(obj.idx_notnan,1)];
%             else
%                 blocks_krig=0:block_size:size(obj.idx_notnan,1);
%                 if blocks_krig(end)~=size(obj.idx_notnan,1)
%                     blocks_krig=[blocks_krig size(obj.idx_notnan,1)];
%                 end
%             end
%             
%             if replications<1
%                 error('Replications must be > 0');
%             end
%             
%             T=obj.stem_model.T;
%             if savetofile
%                 for t=1:T
%                     %create the destination directories
%                     mkdir([directory,'t',num2str(t)]);
%                 end
%                 st_bootkrig_result=[];
%             else
%                 st_bootkrig_result=stem_bootkrig_result(variable_name,obj.stem_grid,obj.stem_model.stem_data.shape);
%             end
%             
%             
%             for b=1:replications
%                 cdisp('Kriging replication ',num2str(b),' started...');
%                 if (b==1)&&not(savetofile)
%                     st_bootkrig_result.Y_hat=nan(obj.Grid_size(1),obj.Grid_size(2),T,rep);
%                 end
%                 
%                 st_par_original=obj.stem_model.stem_par;            
%                 %generate a new st_par randomly
%                 disp('Generating st_par randomly');
%                 st_par_rnd=st_par_original.random(obj.stem_model.stem_EM_result.Hessian);
%                 disp('Generation completed');
%                 obj.stem_model.stem_par=st_par_rnd;
%                 
%                 %do the kriging
%                 st_krig_result = obj.kriging(block_size,variable_name,X,sub,mask,back_transform,no_varcov);
%                
%                 %set the original parameter value
%                 obj.stem_model.stem_par=st_par_original;
%                  
%                 
%                 if not(savetofile)
%                     st_bootkrig_result.Y_hat(:,:,:,b)=st_krig_result.Y_hat;
%                     if not(no_varcov)
%                         st_bootkrig_result.Var_Y_hat(:,:,:,b)=st_krig_result.Var_Y_hat;
%                     end
%                 else
%                     disp('Saving kriging results on HD');
%                     for t=1:T
%                         boot_krig_fixtime.Y_hat=squeeze(st_krig_result.Y_hat(:,:,t));
%                         if not(no_varcov)
%                             boot_krig_fixtime.Var_Y_hat=squeeze(st_krig_result.Var_Y_hat(:,:,t));
%                         end
%                         boot_krig_fixtime.variable_name=variable_name;
%                         boot_krig_fixtime.temporal_index=t;
%                         boot_krig_fixtime.bootstrap_replica_index=b;
%                         boot_krig_fixtime.readme='The actual mask and grid are contained in the file related with t=1 and the first boostrap replica';
%                         if (t==1)&&(b==1)
%                             boot_krig_fixtime.stem_grid=obj.stem_grid;
%                             boot_krig_fixtime.mask=mask;
%                             boot_krig_fixtime.Y=obj.stem_model.stem_data.Y;
%                         else
%                            boot_krig_fixtime.stem_grid=[];
%                            boot_krig_fixtime.mask=[];
%                            boot_krig_fixtime.Y=[];
%                         end
%                         code=num2str(b);
%                         while length(code)<5
%                             code=['0',code];
%                         end
%                         save([directory,'t',num2str(t),'/boot_krig_fixtime_',code],'boot_krig_fixtime');
%                     end
%                     disp('Results saved.');
%                 end
%             end
%          end
        
%         function st_bootkrig_result = bootstrap_kriging(obj,block_krig_size,boot_method,rep,variable_name,X,sub,mask,directory)
%             if nargin<5
%                 error('Not enough input parameters');
%             end
%             if nargin<6
%                 X=[];
%             end
%             if nargin<7
%                 sub=1;
%             end
%             if isempty(sub)==0
%                 sub=1;
%             end
%             if nargin<8
%                 mask=[];
%             else
%                 temp=mask(:);
%                 temp=isnotnan(temp);
%                 obj.idx_notnan=find(isnotnan(temp));
%             end
%             if nargin>8
%                 savetofile=1;
%                 if not(strcmp(directory(end),'/'))
%                     directory=[directory,'/'];
%                 end
%             else
%                 savetofile=0;
%             end
%             
%             if (block_krig_size<0)
%                 error('block_size must be >=0');
%             end
%             if (boot_method<0)||(boot_method>1)
%                 error('boot_method must be either zero or one');
%             end
%             if rep<1
%                 error('rep must be greater than one');
%             end
%             index=obj.stem_model.stem_data.stem_varset.get_index(variable_name);
%             if isempty(index)
%                 error('The variable name is incorrect');
%             end
%             if sub<=0
%                 error('sub must be greater than zero');
%             end
%             
%             Y_original=obj.stem_model.stem_data.Y;
%             residuals=obj.stem_model.stem_EM_result.Y_hat-obj.stem_model.stem_data.Y;
%             %residuals are considered with respect to each site if T>1
%             T=obj.stem_model.T;
%             if savetofile
%                 for t=1:T
%                     %create the destination directories
%                     mkdir([directory,'t',num2str(t)]);
%                 end
%                 st_bootkrig_result=[];
%             else
%                 st_bootkrig_result=stem_bootkrig_result(variable_name,obj.stem_grid,obj.stem_model.stem_data.shape);
%             end
%             if (boot_method==0)&&(T==1)
%                 warning('Not enough temporal steps for site by site bootstrap. Switching to global bootstrap');
%                 boot_method=1;
%             end
%             disp('Residual evaluation...');
%             if boot_method==0
%                 for i=1:size(residuals,1)
%                     L=isnotnan(residuals(i,:));
%                     if T<20
%                         %classic bootstrap
%                         res{i}=residuals(i,L);
%                     else
%                         %block-bootstrap with trend estimation
%                         [a,b]=obj.ma(residuals(i,:),10);
%                         trend{i}=a; %the trend is not restricted to the non-NaN data because it is in a 1:1 relationship with time
%                         res{i}=b(L);
%                     end
%                 end
%             else
%                 L=isnotnan(residuals);
%                 res=residuals(L);
%                 res=res(:);
%             end
%             disp('Residual evaluation ended.');
%             estimated_par=obj.stem_model.stem_par;
%             for b=1:rep
%                 disp(['Bootstrap iteration ',num2str(b),'/',num2str(rep)]);
%                 disp('New data generation...');
%                 if boot_method==0
%                     for i=1:size(obj.stem_model.N,1)
%                         if T<20
%                             %classic bootstrap
%                             indices=round(unifrnd(1,length(res{i}),T,1));
%                             e=res{i}(indices);
%                             %nan are preserved
%                             row=Y_original(i,:)+e;
%                             obj.stem_model.stem_data.set_Y_row(i,row);
%                         else
%                             %block-bootstrap with trend
%                             lres=length(res{i});
%                             if lres<20
%                                 %if the number of actual residuals is less
%                                 %than 20 the block size is defined
%                                 %accordingly
%                                 block_size=floor(lres/2);
%                             else
%                                 %otherwise block_size is 10
%                                 block_size=10;
%                             end
%                             e=[];
%                             while length(e)<T
%                                 %build the residual vector
%                                 idx=round(unifrnd(1,length(res{i})-block_size+1,1,1));
%                                 e=[e res{i}(idx:idx+block_size-1)];
%                             end
%                             if length(e)>T
%                                 %cut the residual vector if it is too long
%                                 e=e(1:T);
%                             end
%                             row=Y_original(i,:)+trend{i}+e; %note the trend
%                             obj.stem_model.stem_data.set_Y_row(i,row);
%                         end
%                     end
%                 else
%                     indices=round(unifrnd(1,length(res),size(obj.stem_model.stem_data.Y,1)*size(obj.stem_model.stem_data.Y,2)));
%                     e=res(indices);
%                     e=reshape(e,size(obj.stem_model.stem_data.Y,1),size(obj.stem_model.stem_data.Y,2));
%                     obj.stem_model.stem_data.set_Y(Y_original+e);
%                     %obj.stem_model.stem_data.Y=obj.stem_model.stem_data.Y+e;
%                 end
%                 %estimation start here
%                 disp('EM estimation');
%                 obj.stem_model.stem_par_initial=estimated_par;
%                 obj.stem_model.EM_estimate(0.001,300,'single');
%                 st_krig_result = obj.kriging(block_krig_size,variable_name,X,sub,mask);
%                 if (b==1)&&not(savetofile)
%                     st_bootkrig_result.Y_hat=nan(obj.Grid_size(1),obj.Grid_size(2),T,rep);
%                 end
%                 %save results
%                 if not(savetofile)
%                     st_bootkrig_result.Y_hat(:,:,:,b)=st_krig_result.Y_hat;
%                     st_bootkrig_result.Var_Y_hat(:,:,:,b)=st_krig_result.Var_Y_hat;
%                 else
%                     for t=1:T
%                         boot_krig_fixtime.Y_hat=squeeze(st_krig_result.Y_hat(:,:,t));
%                         boot_krig_fixtime.Var_Y_hat=squeeze(st_krig_result.Var_Y_hat(:,:,t));
%                         boot_krig_fixtime.Grid_size=obj.Grid_size;
%                         boot_krig_fixtime.variable_name=variable_name;
%                         boot_krig_fixtime.temporal_index=t;
%                         boot_krig_fixtime.bootstrap_replica_index=b;
%                         boot_krig_fixtime.readme='The actual mask and grid are contained in the file related with t=1 and the first boostrap replica';
%                         if (t==1)&&(b==1)
%                             boot_krig_fixtime.Grid=obj.Grid;
%                             boot_krig_fixtime.mask=mask;
%                             boot_krig_fixtime.Y=obj.stem_model.stem_data.Y;
%                         end
%                         code=num2str(b);
%                         while length(code)<5
%                             code=['0',code];
%                         end
%                         save([directory,'t',num2str(t),'/boot_krig_fixtime_',code],'boot_krig_fixtime');
%                     end
%                 end
%                 %reset original data
%                 obj.stem_model.stem_data.update_Y();
%                 delete(st_krig_result);
%             end
%         end
        
        function set.stem_model(obj,stem_model)
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('stem_model must be of class stem_model');
            end
        end
        
        
        function set.X_all(obj,X_all)
            if size(obj.idx_notnan,1)~=size(X_all,1)
                error('The number of rows of X must be equal to the number of non-masked pixel');
            end
            if sum(isnan(X_all(:)))>0
                error('No NaN are allowed in covariates');
            end
            obj.X_all=X_all;
        end
    end
end