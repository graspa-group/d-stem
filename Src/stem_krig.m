%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%% Code website: https://github.com/graspa-group/d-stem                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_krig < handle
    
    %CONSTANTS
    %NN  - number of kriging sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %T   - number of temporal steps
    
    properties
        stem_model=[];          %[stem_model object]        (1x1) the stem_model object of an estimated model
        stem_krig_data=[];      %[stem_krig_data object]    (1x1) the stem_krig_data object
    end
    
    methods
        
        function obj = stem_krig(stem_model,stem_krig_data)
            %DESCRIPTION: constructor of the class stem_krig
            %
            %INPUT
            %
            %stem_model      - [stem_model object]      (1x1) the stem_model object
            %stem_krig_data  - [stem_krig_data object]  (1x1) the stem_krig_data object
            %
            %OUTPUT
            %obj             - [stem_krig object]       (1x1) the stem_krig object      
            
            if nargin<2
                error('All the input arguments must be provided');
            end
            if not(isempty(stem_krig_data.X))
                if not(stem_model.T==size(stem_krig_data.X,3))
                    error('The number of time steps in the matrix X of stem_krig_data must be equal to the number of time steps in the stem_model object');
                end
                
                q=stem_model.stem_par.q;

                for i=1:q
                    if not(isempty(stem_model.stem_data.stem_varset_p.X_bp_name))
                        X_bp_name=stem_model.stem_data.stem_varset_p.X_bp_name{i};
                        if not(isempty(X_bp_name))
                            for j=1:length(X_bp_name)
                                cmp=strcmp(stem_krig_data.X_names,X_bp_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names in stem_krig_data object must include the loading coefficient ',X_bp_name{j}]);
                                end
                            end
                        end
                    end
                    
                    if not(isempty(stem_model.stem_data.stem_varset_p.X_p_name))
                        X_p_name=stem_model.stem_data.stem_varset_p.X_p_name{i};
                        if not(isempty(X_p_name))
                            for j=1:length(X_p_name)
                                cmp=strcmp(stem_krig_data.X_names,X_p_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names in stem_krig_data object must include the loading coefficient ',X_p_name{j}]);
                                end
                            end
                        end
                    end
                    
                    if not(isempty(stem_model.stem_data.stem_varset_p.X_beta_name))
                        X_beta_name=stem_model.stem_data.stem_varset_p.X_beta_name{i};
                        if not(isempty(X_beta_name))
                            for j=1:length(X_beta_name)
                                cmp=strcmp(stem_krig_data.X_names,X_beta_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names in stem_krig_data object must include the loading coefficient ',X_beta_name{j}]);
                                end
                            end
                        end
                    end
                    
                    if not(isempty(stem_model.stem_data.stem_varset_p.X_z_name))
                        X_z_name=stem_model.stem_data.stem_varset_p.X_z_name{i};
                        if not(isempty(X_z_name))
                            for j=1:length(X_z_name)
                                cmp=strcmp(stem_krig_data.X_names,X_z_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names in stem_krig_data object must include the loading coefficient ',X_z_name{j}]);
                                end
                            end
                        end
                    end
                end
            end            
            obj.stem_model=stem_model;
            obj.stem_krig_data=stem_krig_data;
        end
        
        function st_krig_result = kriging(obj,stem_krig_options)
            %Yaqiong
            %{
            if strcmp(stem_krig_options.type,'y-xbeta')&&not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))
                error('The kriging type ''y-xbeta'' can only be used with f-HDGM models');
            end
            %}
            if stem_krig_options.workers==1
                if not(obj.stem_model.stem_par.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                    st_krig_result = obj.kriging_core(stem_krig_options);
                else
                    st_krig_result = obj.kriging_spline_coeff_core(stem_krig_options);
                    %st_krig_result.beta = obj.stem_model.stem_par.beta;
                    %st_krig_result.fda = obj.stem_model.stem_data.stem_fda;
                    %Yaqiong
                    st_krig_result.stem_par = obj.stem_model.stem_par;
                    st_krig_result.stem_fda = obj.stem_model.stem_data.stem_fda;
                end
                %{
                if strcmp(stem_krig_options.type,'y')
                    st_krig_result = obj.kriging_core(stem_krig_options);
                else
                    st_krig_result = obj.kriging_spline_coeff_core(stem_krig_options);
                end
                %}  
            else
                poolobj = parpool(stem_krig_options.workers);
                
                disp('Creating folder for kriging results...');
                folder_name=['kriging_',num2str(round(unifrnd(10000,99999,1,1)))];
                mkdir(folder_name);
                disp(['Temporary kriging results stored in folder ',folder_name]);

                block_krig_size=stem_krig_options.block_size;
                
                stem_krig_options_block=stem_krig_options;
                stem_krig_options_block.block_size=0;
                stem_krig_options_block.workers=1;

                krig_coordinates=obj.stem_krig_data.grid.coordinate;
                
                blocks=0:block_krig_size:size(krig_coordinates,1);
                if not(blocks(end)==size(krig_coordinates,1))
                    blocks=[blocks size(krig_coordinates,1)];
                end
                
                krig_coordinates_block=cell(length(blocks)-1,1);
                for k=1:length(blocks)-1
                    krig_coordinates_block{k}=krig_coordinates((k-1)*block_krig_size+1:k*block_krig_size,:);
                end
                
                disp(['Starting kriging using ',num2str(stem_krig_options.workers),' workers...']);
                %kriging_type=stem_krig_options.type;
                kriging_type = obj.stem_model.stem_par.stem_modeltype;
                parfor k=1:length(blocks)-1
                    disp(['Kriging block: ',num2str(k),' of ',num2str(length(blocks)-1)]);
                    mat_file=matfile([folder_name,'/block',num2str(k,'%010.0f')],'writable',true)

                    obj_stem_krig_grid = stem_grid(krig_coordinates_block{k}, 'deg', 'sparse','point');
                    
                    %Yaqiong
                    if not(strcmp(kriging_type,'f-HDGM'))
                        obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid,X_krig,X_krig_names,mask);
                        obj_stem_krig = stem_krig(obj.stem_model,obj_stem_krig_data);
                        obj_stem_krig_result = obj_stem_krig.kriging_core(stem_krig_options_block);
                    else
                        X_krig=[];
                        X_krig_names=[];
                        mask=[];
                        obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid,X_krig,X_krig_names,mask);
                        obj_stem_krig = stem_krig(obj.stem_model,obj_stem_krig_data);
                        obj_stem_krig_result = obj_stem_krig.kriging_spline_coeff_core(stem_krig_options_block);
                    end
                    
                    mat_file.obj_stem_krig_result_block=obj_stem_krig_result; %save
                end
               
                delete(poolobj);
                
                disp('Merging kriging results...');
                st_krig_result=stem_krig_result(obj.stem_model.stem_data.stem_varset_p.Y_name{1},obj.stem_krig_data.grid,obj.stem_model.stem_data.stem_gridlist_p.grid{1},obj.stem_model.stem_data.shape); 
                st_krig_result.zk_s=[];
                st_krig_result.diag_Pk_s=[];
                st_krig_result.stem_datestamp=obj.stem_model.stem_data.stem_datestamp;
                %Yaqiong
                %st_krig_result.isSplineCoeff=1;
                st_krig_result.stem_par = obj.stem_model.stem_par;
                st_krig_result.stem_fda = obj.stem_model.stem_data.stem_fda;
                
                f=dir([folder_name,'/*.mat']);
                for i=1:length(f)
                    load([folder_name,'/',f(i).name],'obj_stem_krig_result_block');
                    st_krig_result.zk_s=cat(1,st_krig_result.zk_s,obj_stem_krig_result_block.zk_s);
                    st_krig_result.diag_Pk_s=cat(1,st_krig_result.diag_Pk_s,obj_stem_krig_result_block.diag_Pk_s);
                end
                disp(['Removing folder ',folder_name,' with temporary kriging results...']);
                rmdir(folder_name, 's')
                
                if strcmp(obj.stem_krig_data.grid.grid_type,'regular')
                    disp('Kriging result reshaping...');
                    st_krig_result.zk_s=reshape(st_krig_result.zk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result.zk_s,3)]);
                    st_krig_result.diag_Pk_s=reshape(st_krig_result.diag_Pk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result.diag_Pk_s,3)]);
                    disp('Kriging result reshaping ended.');
                    
                    if not(isempty(obj.stem_krig_data.mask))
                        disp('Kriging result masking...');
                        mask=reshape(obj.stem_krig_data.mask,obj.stem_krig_data.grid.grid_size);
                        for t=1:size(st_krig_result.y_hat,3)
                            for h=1:size(st_krig_result.zk_s,3)
                                st_krig_result.zk_s(:,:,t,h)=st_krig_result.zk_s(:,:,t,h).*mask;
                                st_krig_result.diag_Pk_s(:,:,t,h)=st_krig_result.diag_Pk_s(:,:,t,h).*mask;
                            end
                        end
                        disp('Kriging result masking ended.');
                    end
                end
            end
        end
        
        function st_krig_result = kriging_core(obj,stem_krig_options)
            %DESCRIPTION: kriging implementation
            %
            %INPUT
            %
            %obj                - [stem_krig object]            (1x1) the stem_krig object 
            %stem_krig_options  - [stem_krig_options object]    (1x1) a stem_krig_options object

            %OUTPUT
            %st_krig_result     - [stem_krig_result object]     (1x1)              
            
            if nargin<2
                error('Not enough input arguments');
            end

            if not(isa(stem_krig_options,'stem_krig_options'))
                error('The stem_krig_options input argument must be of class stem_krig_options');
            end
            
            if (stem_krig_options.crossval==0)&&isempty(obj.stem_krig_data.X)
                error('X in stem_krig_data cannot be empty');
            end

            if isempty(obj.stem_krig_data.mask)
                idx_notnan=1:size(obj.stem_krig_data.grid.coordinate,1);
                idx_notnan=idx_notnan';
            else
                idx_notnan=find(not(isnan(obj.stem_krig_data.mask)));
            end

            if (stem_krig_options.crossval==1)&&not(obj.stem_model.cross_validation)
                error('The stem_model object does not contain cross-validation information');
            end
            
            if (stem_krig_options.crossval==1)&&not(isempty(obj.stem_krig_data.X))
                error('X in stem_krig_data object must be empty when kriging is performed for cross-validation');
            end
            
            if stem_krig_options.block_size==0
                blocks_krig=[0 size(idx_notnan,1)];
            else
                blocks_krig=0:stem_krig_options.block_size:size(idx_notnan,1);
                if not(blocks_krig(end)==size(idx_notnan,1))
                    blocks_krig=[blocks_krig size(idx_notnan,1)];
                end
            end
            
            K=obj.stem_model.stem_par.k;
            q=obj.stem_model.stem_par.q;
            p=obj.stem_model.stem_par.p;

            %check if X has all the needed covariates
            if stem_krig_options.crossval==0
                idx_bp=cell(q,1);
                idx_p=cell(q,1);
                idx_beta=cell(q,1);
                idx_z=cell(q,1);
                for i=1:q
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_bp_name))
                        X_bp_name=obj.stem_model.stem_data.stem_varset_p.X_bp_name{i};
                        if not(isempty(X_bp_name))
                            idx_bp{i}=[];
                            for j=1:length(X_bp_name)
                                cmp=strcmp(obj.stem_krig_data.X_names,X_bp_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names do not include the covariate ',X_bp_name{j},' included in X_bp']);
                                else
                                    idx_bp{i}=cat(2,idx_bp{i},find(cmp,1));
                                end
                            end
                        else
                            idx_bp{i}=[];
                        end
                    else
                        idx_bp{i}=[];
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_p_name))
                        X_p_name=obj.stem_model.stem_data.stem_varset_p.X_p_name{i};
                        if not(isempty(X_p_name))
                            idx_p{i}=[];
                            for j=1:length(X_p_name)
                                cmp=strcmp(obj.stem_krig_data.X_names,X_p_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names do not include the covariate ',X_p_name{j},' included in X_p']);
                                else
                                    idx_p{i}=cat(2,idx_p{i},find(cmp,1));
                                end
                            end
                        else
                            idx_p{i}=[];
                        end
                    else
                        idx_p{i}=[];
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta_name))
                        X_beta_name=obj.stem_model.stem_data.stem_varset_p.X_beta_name{i};
                        if not(isempty(X_beta_name))
                            idx_beta{i}=[];
                            for j=1:length(X_beta_name)
                                cmp=strcmp(obj.stem_krig_data.X_names,X_beta_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names do not include the covariate ',X_beta_name{j},' included in X_beta']);
                                else
                                    idx_beta{i}=cat(2,idx_beta{i},find(cmp,1));
                                end
                            end
                        else
                            idx_beta{i}=[];
                        end
                    else
                        idx_beta{i}=[];
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z_name))
                        X_z_name=obj.stem_model.stem_data.stem_varset_p.X_z_name{i};
                        if not(isempty(X_z_name))
                            idx_z{i}=[];
                            for j=1:length(X_z_name)
                                cmp=strcmp(obj.stem_krig_data.X_names,X_z_name{j});
                                if sum(cmp)==0
                                    error(['X and X_names do not include the covariate ',X_z_name{j},' included in X_z']);
                                else
                                    idx_z{i}=cat(2,idx_z{i},find(cmp,1));
                                end
                            end
                        else
                            idx_z{i}=[];
                        end
                    else
                        idx_z{i}=[];
                    end
                end
            else
                %kriging using cross-validation data, no need to check X
            end
            
            disp('Kriging started...');

            st_krig_result=cell(q,1);
            for j=1:q
                st_krig_result{j}=stem_krig_result(obj.stem_model.stem_data.stem_varset_p.Y_name{j},obj.stem_krig_data.grid,obj.stem_model.stem_data.stem_gridlist_p.grid{j},obj.stem_model.stem_data.shape);

                if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                    st_krig_result{j}.y_hat=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T);
                    if not(stem_krig_options.no_varcov)
                        st_krig_result{j}.diag_Var_y_hat=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T);
                    end
                    if K>0
                        st_krig_result{j}.E_wp_y1=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T,K);
                    end
                end
                if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                    st_krig_result{j}.zk_s=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T,p);
                    st_krig_result{j}.diag_Pk_s=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T,p);
                end
                st_krig_result{j}.coord_output_block=cell(length(blocks_krig),1);
                st_krig_result{j}.coord_cond_block=cell(length(blocks_krig),1);
            end
            
            for i=1:length(blocks_krig)-1
                ct1=clock;
                if length(blocks_krig)-1>1
                    disp(['Kriging block ',num2str(i),' of ',num2str(length(blocks_krig)-1)]);
                end
                block_krig=(blocks_krig(i)+1):blocks_krig(i+1);
                block_krig_length=length(block_krig);
                
                idx_all=cell(q,1);
                idx_keep=cell(q,1);
                idx_remove=cell(q,1);
                Y_removed=cell(q,1);
                Y_kept=cell(q,1);
                X_bp_removed=cell(q,1);
                X_bp_kept=cell(q,1);
                X_p_removed=cell(q,1);
                X_p_kept=cell(q,1);
                X_beta_removed=cell(q,1);
                X_beta_kept=cell(q,1);
                X_f_removed=cell(q,1); %Yaqiong
                X_f_kept=cell(q,1);
                X_z_removed=cell(q,1);
                X_z_kept=cell(q,1);
                coordinate_removed=cell(q,1);
                coordinate_kept=cell(q,1);
                block_kept_size=zeros(q,1);
                
                isotopic=1;
                for z=2:q
                    if not(size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1)==size(obj.stem_model.stem_data.stem_gridlist_p.grid{z-1}.coordinate,1))
                        isotopic=0;
                    else
                        diff_coord=obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate-obj.stem_model.stem_data.stem_gridlist_p.grid{z-1}.coordinate;
                        if sum(abs(diff_coord(:)))>0
                            isotopic=0;
                        end
                    end
                end
                
                for z=1:q
                    block_coordinates=obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:);
                    
                    if z==1||isotopic==0
                        if stem_krig_options.nn_size>0
                            if strcmp(obj.stem_krig_data.grid.unit,'deg')
                                nn_distance_matrix=zeros(size(obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:),1),size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1));
                                for h=1:size(obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:),1)
                                    nn_distance_matrix(h,:)=distance(block_coordinates(h,1),block_coordinates(h,2),obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(:,1),obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(:,2));
                                end
                            else
                                nn_distance_matrix=pdist2(block_coordinates,obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate);
                            end
                            idx_all{z}=1:size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1);
                            idx_keep{z}=[];
                            for h=1:size(nn_distance_matrix,1)
                                [~,idx]=sort(nn_distance_matrix(h,:));
                                idx=idx(1:stem_krig_options.nn_size);
                                idx_keep{z}=union(idx_keep{z},idx);
                            end
                            idx_keep{z}=idx_keep{z}';
                            idx_remove{z}=setdiff(idx_all{z},idx_keep{z});
                        else
                            idx_all{z}=1:size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1);
                            idx_keep{z}=idx_all{z};
                            idx_remove{z}=[];
                        end
                    else
                        idx_all{z}=idx_all{z-1};
                        idx_keep{z}=idx_keep{z-1};
                        idx_remove{z}=idx_remove{z-1};
                    end
                    block_kept_size(z)=length(idx_keep{z});
                    
                    %Y manage
                    Y_add=nan(block_krig_length,obj.stem_model.T);
                    Y_removed{z}=obj.stem_model.stem_data.stem_varset_p.Y{z}(idx_remove{z},:);
                    Y_kept{z}=obj.stem_model.stem_data.stem_varset_p.Y{z}(idx_keep{z},:);
                    obj.stem_model.stem_data.stem_varset_p.Y{z}=cat(1,Y_kept{z},Y_add);
                    
                    %X manage
                    if stem_krig_options.crossval==0
                        if not(isempty(idx_bp{z}))
                            X_krig_block=obj.stem_krig_data.X(block_krig,idx_bp{z},:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_bp_means{z}(j))/obj.stem_model.stem_data.stem_varset_p.X_bp_stds{z}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_bp{z},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            X_bp_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_bp{z}(idx_remove{z},:,:);
                            X_bp_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_bp{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_bp{z}=cat(1,X_bp_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(idx_p{z}))
                            X_krig_block=obj.stem_krig_data.X(block_krig,idx_p{z},:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_p_means{z}(j))/obj.stem_model.stem_data.stem_varset_p.X_p_stds{z}(j);
                                end
                            end
                            temp=[];
                            for k=1:size(X_krig_block,2)
                                temp1=X_krig_block(:,k,:);
                                temp=cat(4,temp,temp1);
                            end
                            X_krig_block=temp;
                            if size(obj.stem_model.stem_data.stem_varset_p.X_p{z},3)==1
                                X_krig_block=X_krig_block(:,:,1,:);
                            end
                            X_p_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_p{z}(idx_remove{z},:,:);
                            X_p_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_p{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_p{z}=cat(1,X_p_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(idx_beta{z}))
                            X_krig_block=obj.stem_krig_data.X(block_krig,idx_beta{z},:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_beta_means{z}(j))/obj.stem_model.stem_data.stem_varset_p.X_beta_stds{z}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_beta{z},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            X_beta_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_remove{z},:,:);
                            X_beta_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_beta{z}=cat(1,X_beta_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(idx_z{z}))
                            X_krig_block=obj.stem_krig_data.X(block_krig,idx_z{z},:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_z_means{z}(j))/obj.stem_model.stem_data.stem_varset_p.X_z_stds{z}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_z{z},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            X_z_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_remove{z},:,:);
                            X_z_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_z{z}=cat(1,X_z_kept{z},X_krig_block);
                        end
                    else
                        %cross-validation data
                        if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_bp))
                            X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_bp{1}(idx_notnan(block_krig),:,:);
                            X_bp_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_bp{z}(idx_remove{z},:,:);
                            X_bp_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_bp{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_bp{z}=cat(1,X_bp_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_p))
                            X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_p{1}(idx_notnan(block_krig),:,:,:);
                            X_p_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_p{z}(idx_remove{z},:,:);
                            X_p_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_p{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_p{z}=cat(1,X_p_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_beta))
                            X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_beta{1}(idx_notnan(block_krig),:,:);
                            X_beta_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_remove{z},:,:);
                            X_beta_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_beta{z}=cat(1,X_beta_kept{z},X_krig_block);
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_f))
                            X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_f{1}(idx_notnan(block_krig),:,:);
                            X_f_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_f{z}(idx_remove{z},:,:);
                            X_f_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_f{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_f{z}=cat(1,X_f_kept{z},X_krig_block);
                        end %Yaqiong
                        
                        if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_z))
                            X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{z}.X_z{1}(idx_notnan(block_krig),:,:);
                            X_z_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_remove{z},:,:);
                            X_z_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_keep{z},:,:);
                            obj.stem_model.stem_data.stem_varset_p.X_z{z}=cat(1,X_z_kept{z},X_krig_block);
                        end
                    end
                    
                    %Grid manage
                    coordinate_removed{z}=obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(idx_remove{z},:);
                    coordinate_kept{z}=obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(idx_keep{z},:);

                    obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate=cat(1,coordinate_kept{z},block_coordinates);
                    
                    st_krig_result{z}.coord_output_block{i}=block_coordinates;
                    st_krig_result{z}.coord_cond_block{i}=coordinate_kept{z};
                end
                clear Y_add
                clear X_krig_block

                obj.stem_model.stem_data.update_data();
                obj.stem_model.set_distance_matrix('point');
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end
                
                blocks=cumsum(block_kept_size+block_krig_length);
                
                %kriging
                [y_hat,diag_Var_y_hat,E_wp_y1,diag_Var_wp_y1,stem_kalmansmoother_result]=obj.E_step_kriging(stem_krig_options.no_varcov,stem_krig_options.crossval);
                for j=1:q
                    if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                        st_krig_result{j}.y_hat(idx_notnan(block_krig),:)=y_hat(blocks(j)-block_krig_length+1:blocks(j),:);
                        if not(stem_krig_options.no_varcov)
                            st_krig_result{j}.diag_Var_y_hat(idx_notnan(block_krig),:)=diag_Var_y_hat(blocks(j)-block_krig_length+1:blocks(j),:);
                            if K>0
                                st_krig_result{j}.diag_Var_wp_y1(idx_notnan(block_krig),:,:)=diag_Var_wp_y1(blocks(j)-block_krig_length+1:blocks(j),:,:);
                            end
                        end
                        if K>0
                            st_krig_result{j}.E_wp_y1(idx_notnan(block_krig),:,:)=E_wp_y1(blocks(j)-block_krig_length+1:blocks(j),:,:);
                        end
                    end
                    
                    if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                        if obj.stem_model.stem_data.stem_modeltype.is('f-HDGM')
                            temp=block_kept_size(1)*ones(1,p);
                            blocks_z=cumsum(temp+block_krig_length);
                        else
                            blocks_z=blocks;
                        end
                        for h=1:size(st_krig_result{j}.zk_s,3)
                            st_krig_result{j}.zk_s(idx_notnan(block_krig),:,h)=stem_kalmansmoother_result.zk_s(blocks_z(h)-block_krig_length+1:blocks_z(h),2:end);
                            for t=1:size(st_krig_result{j}.diag_Pk_s,2)
                                st_krig_result{j}.diag_Pk_s(idx_notnan(block_krig),t,h)=diag(stem_kalmansmoother_result.Pk_s{t+1}(blocks_z(h)-block_krig_length+1:blocks_z(h),blocks_z(h)-block_krig_length+1:blocks_z(h)));
                            end
                        end
                    end
                end
                
                %restore original
                for j=1:q
                    obj.stem_model.stem_data.stem_varset_p.Y{j}(end-block_krig_length+1:end,:)=[];
                    obj.stem_model.stem_data.stem_varset_p.Y{j}([idx_keep{j} idx_remove{j}],:)=[Y_kept{j};Y_removed{j}];
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_bp))
                        obj.stem_model.stem_data.stem_varset_p.X_bp{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_bp{j}([idx_keep{j} idx_remove{j}],:,:)=[X_bp_kept{j};X_bp_removed{j}];
                    end
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_p))
                        obj.stem_model.stem_data.stem_varset_p.X_p{j}(end-block_krig_length+1:end,:,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_p{j}([idx_keep{j} idx_remove{j}],:,:,:)=[X_p_kept{j};X_p_removed{j}];
                    end
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta))
                        obj.stem_model.stem_data.stem_varset_p.X_beta{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_beta{j}([idx_keep{j} idx_remove{j}],:,:)=[X_beta_kept{j};X_beta_removed{j}];
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_f))
                        obj.stem_model.stem_data.stem_varset_p.X_f{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_f{j}([idx_keep{j} idx_remove{j}],:,:)=[X_f_kept{j};X_f_removed{j}];
                    end %Yaqiong
                    
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z))
                        obj.stem_model.stem_data.stem_varset_p.X_z{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_z{j}([idx_keep{j} idx_remove{j}],:,:)=[X_z_kept{j};X_z_removed{j}];
                    end
                    
                    obj.stem_model.stem_data.stem_gridlist_p.grid{j}.coordinate(end-block_krig_length+1:end,:)=[];
                    obj.stem_model.stem_data.stem_gridlist_p.grid{j}.coordinate([idx_keep{j} idx_remove{j}],:)=[coordinate_kept{j};coordinate_removed{j}];
                end
                
                obj.stem_model.stem_data.update_data(); 
                obj.stem_model.set_distance_matrix('point');
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end
                
                ct2=clock;
                disp(['Kriging block ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end

            if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                if stem_krig_options.back_transform&&(obj.stem_model.stem_data.stem_varset_p.standardized||obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    disp('Back-transformation...');
                    for j=1:q
                        s=obj.stem_model.stem_data.stem_varset_p.Y_stds{j};
                        m=obj.stem_model.stem_data.stem_varset_p.Y_means{j};
                        if (obj.stem_model.stem_data.stem_varset_p.standardized)&&not(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                            st_krig_result{j}.y_hat=st_krig_result.y_hat*s+m;
                            if not(stem_krig_options.no_varcov)
                                st_krig_result{j}.diag_Var_y_hat=st_krig_result.diag_Var_y_hat*s^2;
                            end
                        end
                        
                        if (obj.stem_model.stem_data.stem_varset_p.standardized)&&(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                            y_hat=st_krig_result{j}.y_hat;
                            var_y_hat=st_krig_result{j}.diag_Var_y_hat;
                            st_krig_result{j}.y_hat=exp(y_hat*s+m+(var_y_hat*s^2)/2);
                            if not(stem_krig_options.no_varcov)
                                st_krig_result{j}.diag_Var_y_hat=(exp(var_y_hat*s^2)-1).*exp(2*(y_hat*s+m)+(var_y_hat*s^2));
                            end
                        end
                    end
                    disp('Back-transformation ended.');
                end
            end
            
            if strcmp(obj.stem_krig_data.grid.grid_type,'regular')
                disp('Data reshaping...');
                for j=1:q
                    if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                        st_krig_result{j}.y_hat=reshape(st_krig_result{j}.y_hat,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T]);
                        if not(stem_krig_options.no_varcov)
                            st_krig_result{j}.diag_Var_y_hat=reshape(st_krig_result{j}.diag_Var_y_hat,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T]);
                        end
                        if K>0
                            st_krig_result{j}.E_wp_y1=reshape(st_krig_result{j}.E_wp_y1,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,K]);
                        end
                    end
                    if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                        st_krig_result{j}.zk_s=reshape(st_krig_result{j}.zk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result{j}.zk_s,3)]);
                        st_krig_result{j}.diag_Pk_s=reshape(st_krig_result{j}.diag_Pk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result{j}.diag_Pk_s,3)]);
                    end
                end
                disp('Date reshaping ended.');
                
                if not(isempty(obj.stem_krig_data.mask))
                    disp('Applying mask...');
                    mask=reshape(obj.stem_krig_data.mask,obj.stem_krig_data.grid.grid_size);
                    for j=1:q
                        for t=1:size(st_krig_result{j}.y_hat,3)
                            if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))||stem_krig_options.crossval==1
                                st_krig_result{j}.y_hat(:,:,t)=st_krig_result{j}.y_hat(:,:,t).*mask;
                                if not(no_varcov)
                                    st_krig_result{j}.diag_Var_y_hat(:,:,t)=st_krig_result{j}.diag_Var_y_hat(:,:,t).*mask;
                                end
                                for k=1:K
                                    st_krig_result{j}.E_wp_y1(:,:,t,k)=st_krig_result{j}.E_wp_y1(:,:,t,k).*mask;
                                end
                            end
                            if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                                for h=1:size(st_krig_result{j}.zk_s,3)
                                    st_krig_result{j}.zk_s(:,:,t,h)=st_krig_result{j}.zk_s(:,:,t,h).*mask;
                                    st_krig_result{j}.diag_Pk_s(:,:,t,h)=st_krig_result{j}.diag_Pk_s(:,:,t,h).*mask;
                                end
                            end
                        end
                    end
                    disp('Mask applied.');
                end
            end
            
            for j=1:q
                st_krig_result{j}.variable_name=obj.stem_model.stem_data.stem_varset_p.Y_name{j};
                st_krig_result{j}.stem_datestamp=obj.stem_model.stem_data.stem_datestamp;
            end
        end
        
        function st_krig_result = kriging_spline_coeff_core(obj,stem_krig_options)
            %DESCRIPTION: kriging of the spline coefficients. Only works with f-HDGM models
            %
            %INPUT
            %
            %obj                - [stem_krig object]            (1x1) the stem_krig object 
            %stem_krig_options  - [stem_krig_options object]    (1x1) a stem_krig_options object

            %OUTPUT
            %st_krig_result     - [stem_krig_result object]     (1x1)  
            
            if not(obj.stem_model.stem_data.stem_modeltype.is('f-HDGM'))
                error('This method only works with f-HDGM models');
            end
            
            if nargin<2
                error('Not enough input arguments');
            end

            if not(isa(stem_krig_options,'stem_krig_options'))
                error('The stem_krig_options input argument must be of class stem_krig_options');
            end
            
            if isempty(obj.stem_krig_data.mask)
                idx_notnan=1:size(obj.stem_krig_data.grid.coordinate,1);
                idx_notnan=idx_notnan';
            else
                idx_notnan=find(not(isnan(obj.stem_krig_data.mask)));
            end
            
            if stem_krig_options.block_size==0
                blocks_krig=[0 size(idx_notnan,1)];
            else
                blocks_krig=0:stem_krig_options.block_size:size(idx_notnan,1);
                if not(blocks_krig(end)==size(idx_notnan,1))
                    blocks_krig=[blocks_krig size(idx_notnan,1)];
                end
            end
            
            q=obj.stem_model.stem_par.q;
            p=obj.stem_model.stem_par.p;

            disp('Kriging started...');

            %the stem_krig_result object is created and initialized
            st_krig_result=stem_krig_result(obj.stem_model.stem_data.stem_varset_p.Y_name{1},obj.stem_krig_data.grid,obj.stem_model.stem_data.stem_gridlist_p.grid{1},obj.stem_model.stem_data.shape); %note {1} since names and grids are equal when model_type is f-HDGM
            st_krig_result.zk_s=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T,p);
            st_krig_result.diag_Pk_s=zeros(size(obj.stem_krig_data.grid.coordinate,1),obj.stem_model.T,p);
            st_krig_result.coord_output_block=cell(length(blocks_krig),1);
            st_krig_result.coord_cond_block=cell(length(blocks_krig),1);
            st_krig_result.stem_datestamp=obj.stem_model.stem_data.stem_datestamp;
            %st_krig_result.isSplineCoeff=1;
            
            for i=1:length(blocks_krig)-1
                ct1=clock;
                if length(blocks_krig)-1>1
                    disp(['Kriging block ',num2str(i),' of ',num2str(length(blocks_krig)-1)]);
                end
                block_krig=(blocks_krig(i)+1):blocks_krig(i+1);
                block_krig_length=length(block_krig);
                
                idx_all=cell(q,1);
                idx_keep=cell(q,1);
                idx_remove=cell(q,1);
                Y_removed=cell(q,1);
                Y_kept=cell(q,1);
                X_beta_removed=cell(q,1);
                X_beta_kept=cell(q,1);
                X_f_removed=cell(q,1);
                X_f_kept=cell(q,1); %Yaqiong
                X_z_removed=cell(q,1);
                X_z_kept=cell(q,1);
                coordinate_removed=cell(q,1);
                coordinate_kept=cell(q,1);
                block_kept_size=zeros(q,1);

                for z=1:q
                    block_coordinates=obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:);

                    if z==1
                        if stem_krig_options.nn_size>0
                            if strcmp(obj.stem_krig_data.grid.unit,'deg')
                                nn_distance_matrix=zeros(size(obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:),1),size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1));
                                for h=1:size(obj.stem_krig_data.grid.coordinate(idx_notnan(block_krig),:),1)
                                    nn_distance_matrix(h,:)=distance(block_coordinates(h,1),block_coordinates(h,2),obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(:,1),obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(:,2));
                                end
                            else
                                nn_distance_matrix=pdist2(block_coordinates,obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate);
                            end
                            idx_all{z}=1:size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1);
                            idx_keep{z}=[];
                            for h=1:size(nn_distance_matrix,1)
                                [~,idx]=sort(nn_distance_matrix(h,:));
                                idx=idx(1:stem_krig_options.nn_size);
                                idx_keep{z}=union(idx_keep{z},idx);
                            end
                            idx_keep{z}=idx_keep{z}';
                            idx_remove{z}=setdiff(idx_all{z},idx_keep{z});
                        else
                            idx_all{z}=1:size(obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate,1);
                            idx_keep{z}=idx_all{z};
                            idx_remove{z}=[];
                        end
                    else
                        %This can be done since f-HDGM models are isotopic
                        idx_all{z}=idx_all{z-1};
                        idx_keep{z}=idx_keep{z-1};
                        idx_remove{z}=idx_remove{z-1};
                    end
                    block_kept_size(z)=length(idx_keep{z});

                    %Y manage
                    Y_add=nan(block_krig_length,obj.stem_model.T);
                    Y_removed{z}=obj.stem_model.stem_data.stem_varset_p.Y{z}(idx_remove{z},:);
                    Y_kept{z}=obj.stem_model.stem_data.stem_varset_p.Y{z}(idx_keep{z},:);
                    obj.stem_model.stem_data.stem_varset_p.Y{z}=cat(1,Y_kept{z},Y_add);

                    %X manage
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta))
                        if size(obj.stem_model.stem_data.stem_varset_p.X_beta{z},3)==1
                            X_krig_block=zeros(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_beta{z},2));                        
                        else
                            X_krig_block=zeros(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_beta{z},2),size(obj.stem_model.stem_data.stem_varset_p.X_beta{z},3));
                        end
                        X_beta_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_remove{z},:,:);
                        X_beta_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_beta{z}(idx_keep{z},:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_beta{z}=cat(1,X_beta_kept{z},X_krig_block);
                    end
                    %Yaqiong
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_f))
                        if size(obj.stem_model.stem_data.stem_varset_p.X_f{z},3)==1
                            %note that nan instead of zeros
                            X_krig_block=nan(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_f{z},2));                        
                        else
                            X_krig_block=nan(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_f{z},2),size(obj.stem_model.stem_data.stem_varset_p.X_f{z},3));
                        end
                        X_f_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_f{z}(idx_remove{z},:,:);
                        X_f_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_f{z}(idx_keep{z},:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_f{z}=cat(1,X_f_kept{z},X_krig_block);
                    end

                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z))
                        if size(obj.stem_model.stem_data.stem_varset_p.X_z{z},3)==1
                            X_krig_block=zeros(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_z{z},2));                        
                        else
                            X_krig_block=zeros(length(block_krig),size(obj.stem_model.stem_data.stem_varset_p.X_z{z},2),size(obj.stem_model.stem_data.stem_varset_p.X_z{z},3));
                        end
                        X_z_removed{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_remove{z},:,:);
                        X_z_kept{z}=obj.stem_model.stem_data.stem_varset_p.X_z{z}(idx_keep{z},:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_z{z}=cat(1,X_z_kept{z},X_krig_block);
                    end

                    %Grid manage
                    coordinate_removed{z}=obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(idx_remove{z},:);
                    coordinate_kept{z}=obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate(idx_keep{z},:);

                    obj.stem_model.stem_data.stem_gridlist_p.grid{z}.coordinate=cat(1,coordinate_kept{z},block_coordinates);
                    
                    st_krig_result.coord_output_block{i}=block_coordinates;
                    st_krig_result.coord_cond_block{i}=coordinate_kept{z};
                end
                clear Y_add
                clear X_krig_block

                obj.stem_model.stem_data.update_data();
                obj.stem_model.set_distance_matrix('point');
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end
                                
                %kriging
                [~,~,~,~,stem_kalmansmoother_result]=obj.E_step_kriging(stem_krig_options.no_varcov,stem_krig_options.crossval);
                    
                %save results in the stem_krig_result object
                temp=block_kept_size(1)*ones(1,p);
                blocks_z=cumsum(temp+block_krig_length);
                for h=1:size(st_krig_result.zk_s,3)
                    st_krig_result.zk_s(idx_notnan(block_krig),:,h)=stem_kalmansmoother_result.zk_s(blocks_z(h)-block_krig_length+1:blocks_z(h),2:end);
                    for t=1:size(st_krig_result.diag_Pk_s,2)
                        st_krig_result.diag_Pk_s(idx_notnan(block_krig),t,h)=diag(stem_kalmansmoother_result.Pk_s{t+1}(blocks_z(h)-block_krig_length+1:blocks_z(h),blocks_z(h)-block_krig_length+1:blocks_z(h)));
                    end
                end
             
                %restore original
                for j=1:q
                    obj.stem_model.stem_data.stem_varset_p.Y{j}(end-block_krig_length+1:end,:)=[];
                    obj.stem_model.stem_data.stem_varset_p.Y{j}([idx_keep{j} idx_remove{j}],:)=[Y_kept{j};Y_removed{j}];
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta))
                        obj.stem_model.stem_data.stem_varset_p.X_beta{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_beta{j}([idx_keep{j} idx_remove{j}],:,:)=[X_beta_kept{j};X_beta_removed{j}];
                    end
                    %Yaqiong
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_f))
                        obj.stem_model.stem_data.stem_varset_p.X_f{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_f{j}([idx_keep{j} idx_remove{j}],:,:)=[X_f_kept{j};X_f_removed{j}];
                    end
                    if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z))
                        obj.stem_model.stem_data.stem_varset_p.X_z{j}(end-block_krig_length+1:end,:,:)=[];
                        obj.stem_model.stem_data.stem_varset_p.X_z{j}([idx_keep{j} idx_remove{j}],:,:)=[X_z_kept{j};X_z_removed{j}];
                    end
                    
                    obj.stem_model.stem_data.stem_gridlist_p.grid{j}.coordinate(end-block_krig_length+1:end,:)=[];
                    obj.stem_model.stem_data.stem_gridlist_p.grid{j}.coordinate([idx_keep{j} idx_remove{j}],:)=[coordinate_kept{j};coordinate_removed{j}];
                end
                
                obj.stem_model.stem_data.update_data();
                obj.stem_model.set_distance_matrix('point');
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end

                ct2=clock;
                disp(['Kriging block ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end

            if strcmp(obj.stem_krig_data.grid.grid_type,'regular')
                disp('Data reshaping...');
                st_krig_result.zk_s=reshape(st_krig_result.zk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result.zk_s,3)]);
                st_krig_result.diag_Pk_s=reshape(st_krig_result.diag_Pk_s,[obj.stem_krig_data.grid.grid_size,obj.stem_model.T,size(st_krig_result.diag_Pk_s,3)]);
                disp('Date reshaping ended.');
                
                if not(isempty(obj.stem_krig_data.mask))
                    disp('Applying mask...');
                    mask=reshape(obj.stem_krig_data.mask,obj.stem_krig_data.grid.grid_size);
                    for t=1:size(st_krig_result.y_hat,3)
                        for h=1:size(st_krig_result.zk_s,3)
                            st_krig_result.zk_s(:,:,t,h)=st_krig_result.zk_s(:,:,t,h).*mask;
                            st_krig_result.diag_Pk_s(:,:,t,h)=st_krig_result.diag_Pk_s(:,:,t,h).*mask;
                        end
                    end
                    disp('Mask applied.');
                end
            end
        end

        function [y_hat,diag_Var_y_hat,E_wp_y1,diag_Var_wp_y1,st_kalmansmoother_result] = E_step_kriging(obj,no_varcov,crossval)
            %DESCRIPTION: kriging is based on the E-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_krig object]                    (1x1)
            %no_varcov                      - [boolean]                             (1x1) 1: the variance of the kriged variable is not computed; 0: the variance is computed;
            %crossval                       - [boolean]                             (1x1) 1: kriging is done for cross validation; 0: actual kriging
            %
            %OUTPUT
            %y_hat                          - [double]                              (NNxT) the variable estimated over the kriging sites
            %diag_Var_y_hat                 - [double]                              (NNxT) the variance of the variable estimated over the kriging sites
            %E_wp_y1                        - [double]                              (NNxTxK) E[wp|Y(1)] over the kriging sites
            %diag_Var_wp_y1                 - [double]                              (NNxTxK) diagonals of Var[wp|Y(1)] over the kriging sites
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object]   (1x1) the stem_kalmansmoother_result object
            
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            T=obj.stem_model.stem_data.T;
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            if p>0
                if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                    st_kalman=stem_kalman(obj.stem_model);
                    compute_logL=0;
                    enable_varcov_computation=0;
                    block_tapering_block_size=0;
                    %note that, for kriging, block tapering is disabled despite what happened in model estimation
                    [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,M] = st_kalman.smoother(compute_logL,enable_varcov_computation,[],[],block_tapering_block_size);
                else
                    [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,~,~,aj_bp,M] = obj.stem_model.get_sigma();
                    st_kalmansmoother_result=obj.stem_model.stem_EM_result.stem_kalmansmoother_result;
                end
                rr=size(sigma_Z,1);
                if not(obj.stem_model.stem_data.X_tv)
                    if obj.stem_model.stem_data.stem_modeltype.is('HDGM')
                        temp=obj.stem_model.stem_data.X_z{1};
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[obj.stem_model.stem_data.X_z{1};zeros(N-size(obj.stem_model.stem_data.X_z{1},1),size(obj.stem_model.stem_data.X_z{1},2))];
                    end
                    %X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    
                    if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                        if obj.stem_model.tapering
                            %is it possible to improve the sparse matrix var_Zt?
                            var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                        else
                            var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                        end
                    end
                    if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                        var_Yt=sigma_geo+var_Zt;
                    end
                end
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,~,~,aj_bp,M] = obj.stem_model.get_sigma();
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[],[],[]);
                var_Zt=[];
                if not(obj.stem_model.stem_data.X_tv)
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
                rr=0;
            end
            
            if obj.stem_model.stem_data.stem_modeltype.is('f-HDGM')&&crossval==0
                %no need to compute the rest, only st_kalmansmoother_result is needed
                y_hat=[];
                diag_Var_y_hat=[];
                E_wp_y1=[];
                diag_Var_wp_y1=[];
                return
            end

            res=data.Y;
            res(isnan(res))=0;
            if not(isempty(data.X_beta))
                Xbeta=zeros(N,T);
                if data.X_beta_tv
                    for t=1:T
                        if size(data.X_beta{t},1)<N
                            X_beta_orlated=[data.X_beta{t};zeros(N-size(data.X_beta{t},1),size(data.X_beta{t},2))];
                        else
                            X_beta_orlated=data.X_beta{t};
                        end
                        Xbeta(:,t)=X_beta_orlated*par.beta;
                    end
                else
                    if size(data.X_beta{1},1)<N
                        X_beta_orlated=[data.X_beta{1};zeros(N-size(data.X_beta{1},1),size(data.X_beta{1},2))];
                    else
                        X_beta_orlated=data.X_beta{1};
                    end
                    Xbeta=repmat(X_beta_orlated*par.beta,1,T);
                end
                res=res-Xbeta;
                y_hat=Xbeta;
            else
                y_hat=zeros(size(res));
            end

            if not(no_varcov)
                diag_Var_e_y1=zeros(N,T);
            else
                diag_Var_e_y1=[];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if not(isempty(data.X_bp))
                %cov_wb_yz time invariant case
                if not(data.X_bp_tv)
                    cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),data.X_bp{1},'r'),aj_bp,'r');
                end
                E_wb_y1=zeros(Nb,T);
                diag_Var_wb_y1=zeros(Nb,T);
                cov_wb_z_y1=zeros(Nb,rr,T);
            end
            
            if not(isempty(data.X_p))
                aj_p=[ones(N-Nb,1); zeros(Nb,1)];
            end

            if not(isempty(data.X_p))
                if not(data.X_p_tv)
                    cov_wp_y=cell(K,1);
                    %cov_wp_yz time invariant case
                    for k=1:K
                        cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p{1}(:,k),'r'),aj_p,'r');
                        %cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p{1}(:,k),'r'),aj_p(:,k),'r');
                    end
                end
                cov_wpk_wph_y1=cell(K,K);
                for h=1:K
                    for k=h+1:K
                        cov_wpk_wph_y1{k,h}=zeros(Np,T);
                    end
                end
                E_wp_y1=zeros(Np,T,K);
                if not(no_varcov)
                    diag_Var_wp_y1=zeros(Np,T,K);
                else
                    diag_Var_wp_y1=[];
                end
                cov_wp_z_y1=zeros(Np,rr,T,K);
            else
                E_wp_y1=[];
                diag_Var_wp_y1=[];
            end

            if not(isempty(data.X_bp)) && not(isempty(data.X_p))
                M_cov_wb_wp_y1=zeros(N,T,K);
            else
                M_cov_wb_wp_y1=[];
            end

            for t=1:T
                %missing at time t
                Lt=not(isnan(data.Y(:,t)));

                if data.X_bp_tv
                    tBP=t;
                else
                    tBP=1;
                end
                if data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_p_tv
                    tP=t;
                else
                    tP=1;
                end

                %evaluate var_yt in the time variant case
                if data.X_tv
                    if not(isempty(obj.stem_model.stem_data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_model.stem_data.X_bp{tBP},'b'),aj_bp,'b');
                    end
                    if not(isempty(obj.stem_model.stem_data.X_p))
                        if isempty(sigma_geo)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                       
                        for k=1:size(obj.stem_model.stem_data.X_p{1},2)
                            %sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p{tP}(:,k),'b'),aj_p(:,k),'b');
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p{tP}(:,k),'b'),aj_p,'b');
                        end
                    end
                    if isempty(sigma_geo)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end

                    if p>0
                        if obj.stem_model.stem_data.stem_modeltype.is({'HDGM','f-HDGM'})
                            temp=obj.stem_model.stem_data.X_z{tT};
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[obj.stem_model.stem_data.X_z{tT};zeros(N-size(obj.stem_model.stem_data.X_z{tT},1),size(obj.stem_model.stem_data.X_z{tT},2))];
                        end
                        %X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');

                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            if obj.stem_model.tapering
                                %is it possible to improve the sparse matrix var_Zt?
                                var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                            else
                                var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                            end
                            if isempty(sigma_geo)
                                var_Yt=var_Zt;
                            else
                                var_Yt=sigma_geo+var_Zt;
                            end
                        else
                            var_Yt=[];
                            var_Zt=[];
                        end
                    else
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            var_Yt=sigma_geo;
                        else
                            var_Yt=[];
                        end
                    end
                end

                %check if the temporal loadings are time variant
                if p>0
                    if obj.stem_model.product_step>0
                        blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(X_z_orlated(blocks(i)+1:blocks(i+1),:)*st_kalmansmoother_result.Pk_s{t+1}*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                        end
                    else
                        temp=X_z_orlated*st_kalmansmoother_result.Pk_s{t+1};
                        diag_Var_e_y1(:,t)=diag(temp*X_z_orlated');
                    end
                    %update E(e|y1)
                    temp=st_kalmansmoother_result.zk_s(:,t+1);
                    y_hat(:,t)=y_hat(:,t)+X_z_orlated*temp;
                end

                if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), X_z_orlated(Lt,:)*sigma_Z; sigma_Z*X_z_orlated(Lt,:)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end

                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=stem_misc.chol_solve(chol_H1t,temp2(r));
                        clear temp2
                    else
                        chol_H1t=chol(H1t);
                        cs=stem_misc.chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end

                if not(isempty(data.X_bp))
                    %check if the pixel loadings are time variant
                    if data.X_bp_tv
                        %cov_wb_yz time variant case
                        cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),data.X_bp(:,1,tBP),'r'),aj_bp,'r');
                    end
                    cov_wb_y1z=[cov_wb_y(:,Lt),zeros(size(cov_wb_y,1),rr)];

                    %compute E(w_b|y1);
                    E_wb_y1(:,t)=cov_wb_y1z*cs;

                    if not(no_varcov)
                        %compute diag(Var(w_b|y1))
                        if obj.stem_model.tapering
                            temp_b(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wb_y1z(:,r)');
                        else
                            temp_b=stem_misc.chol_solve(chol_H1t,cov_wb_y1z');
                        end

                        blocks=0:80:size(diag_Var_wb_y1,1);
                        if not(blocks(end)==size(diag_Var_wb_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_wb_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)=diag(sigma_W_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1))-cov_wb_y1z(blocks(i)+1:blocks(i+1),:)*temp_b(:,blocks(i)+1:blocks(i+1)));
                        end
                    end

                    if (p>0)&&(not(no_varcov))
                        %compute cov(w_b,z|y1)
                        cov_wb_z_y1(:,:,t)=temp_b(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s{t+1};
                        blocks=0:80:size(diag_Var_wb_y1,1);
                        if not(blocks(end)==size(diag_Var_wb_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_wb_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)+diag(cov_wb_z_y1(blocks(i)+1:blocks(i+1),:,t)*temp_b(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                        end
                        clear temp_b
                        %update diag(Var(e|y1))
                        temp=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                        if obj.stem_model.product_step>0
                            blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %notare 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                        end
                    else
                        cov_wb_z_y1=[];
                        clear temp_b
                    end
                    %update y_hat
                    y_hat(:,t)=y_hat(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                    %update diag(Var(e|y1))
                    if not(no_varcov)
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(diag_Var_wb_y1(:,t),M,'l'),data.X_bp(:,1,tBP),'b'),aj_bp,'b'); %tested
                    end
                end

                if not(isempty(data.X_p))
                    %check if the point loadings are time variant
                    if data.X_p_tv
                        %cov_wp_yz time invariant case
                        
                        for k=1:K
                            %cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p(:,1,tP,k),'r'),aj_p(:,k),'r');
                            cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p{tP}(:,k),'r'),aj_p,'r');
                        end
                    end
                    temp_p=cell(K,1);
                    for k=1:K
                        cov_wp_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                        %compute E(w_p_k|y1);
                        E_wp_y1(:,t,k)=cov_wp_y1z*cs;

                        if not(no_varcov)
                            %compute diag(Var(w_p_k|y1))
                            if obj.stem_model.tapering
                                temp_p{k}(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wp_y1z(:,r)');
                            else
                                temp_p{k}=stem_misc.chol_solve(chol_H1t,cov_wp_y1z');
                            end

                            blocks=0:80:size(diag_Var_wp_y1(:,t,k),1);
                            if not(blocks(end)==size(diag_Var_wp_y1(:,t,k),1))
                                blocks=cat(2,blocks,size(diag_Var_wp_y1(:,t,k),1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(sigma_W_p{k}(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1))-cov_wp_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                            end
                        end

                        if (p>0)&&(not(no_varcov))
                            %compute cov(w_p,z|y1)
                            cov_wp_z_y1(:,:,t,k)=temp_p{k}(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s{t+1};
                            blocks=0:80:size(diag_Var_wp_y1(:,t,k),1);
                            if not(blocks(end)==size(diag_Var_wp_y1(:,t,k),1))
                                blocks=cat(2,blocks,size(diag_Var_wp_y1(:,t,k),1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)+diag(cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                            end
                            %update diag(Var(e|y1))
                            
                            %temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                            temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),data.X_p{tP}(:,k),'l'),aj_p,'l');
                            if obj.stem_model.product_step>0
                                blocks=0:obj.stem_model.product_step:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %notare 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                            end
                        else
                            cov_wp_z_y1=[];
                        end
                        %y_hat
                        
                        y_hat(:,t)=y_hat(:,t)+stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,k),data.X_p{tP}(:,k),'l'),aj_p,'l');

                        if not(no_varcov)
                            %update diag(Var(e|y1))
                            
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(diag_Var_wp_y1(:,t,k),data.X_p{tP}(:,k),'b'),aj_p,'b'); %K varianze

                            if not(isempty(data.X_bp))
                                %compute M_cov(w_b,w_p|y1); cio? M*cov(w_b,w_p|y1) da tenere in considerazione nelle forme chiuse!
                                if obj.stem_model.product_step>0
                                    blocks=0:obj.stem_model.product_step:length(M);
                                    if not(blocks(end)==length(M))
                                        blocks=cat(2,blocks,length(M));
                                    end
                                    for i=1:length(blocks)-1
                                        %tested
                                        if p>0
                                            M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1))+cov_wb_z_y1(M(blocks(i)+1:blocks(i+1)),:,t)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1))); %ha gia' l'stem_misc.M_apply su left!!
                                        else
                                            M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if p>0
                                        M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M))+cov_wb_z_y1(M,:,t)*temp_p{k}(end-rr+1:end,1:length(M))); %ha gi? l'stem_misc.M_apply su left!!
                                    else
                                        M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M)));
                                    end
                                end
                                %update diag(Var(e|y1)) - tested
                                temp=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),data.X_bp{tBP}(:,1),'l'),aj_bp,'l');
                                aj_p=[ones(N-Nb,1); zeros(Nb,1)];
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[data.X_p{tP}(:,k);zeros(Nb,1)],'l'),aj_p,'l');
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    end

                    if (K>1)&&(not(no_varcov))
                        %compute cov(w_pk,w_ph|y1);
                        for h=1:K
                            for k=h+1:K
                                cov_wpk_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                                if obj.stem_model.product_step>0
                                    blocks=0:obj.stem_model.product_step:size(cov_wpk_y1z,1);
                                    if not(blocks(end)==size(cov_wpk_y1z,1))
                                        blocks=cat(2,blocks,size(cov_wpk_y1z,1));
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wp_z_y1))
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1))+cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{h}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wp_z_y1))
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h}+cov_wp_z_y1(:,:,t,k)*temp_p{h}(end-rr+1:end,:));
                                    else
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h});
                                    end
                                end
                               
                                %temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                %temp=stem_misc.D_apply(stem_misc.D_apply(temp,[data.X_p(:,1,tP,h);zeros(Nb,1)],'l'),aj_p(:,h),'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),data.X_p{tP}(:,k),'l'),aj_p,'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[data.X_p{tP}(:,h);zeros(Nb,1)],'l'),aj_p,'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    end
                    clear temp_p
                end
                if data.X_tv
                    sigma_geo=[];
                end
            end
            diag_Var_y_hat=diag_Var_e_y1;             
        end
          
        %Class set methods
        function set.stem_model(obj,stem_model)
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('stem_model must be of class stem_model');
            end
        end
        
        function set.stem_krig_data(obj,stem_krig_data)
            if isa(stem_krig_data,'stem_krig_data')
                obj.stem_krig_data=stem_krig_data;
            else
                error('stem_krig_data must be of class stem_krig_data');
            end
        end
        
        
    end
end