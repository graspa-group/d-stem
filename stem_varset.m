%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_varset < handle
    properties
        Y={};               %{q x 1} cells. Each cell c_i [n_i x T double] Variable data
        X_rg={};            %{q x 1} cells. Each cell c_i [n_i x 1 x T|1 double] Loading
        X_beta={};          %{q x 1} cells. Each cell c_i [n_i x n_beta_i x T|1 double]
        X_time={};          %{q x 1} cells. Each cell c_i [n_i x n_time_i x T|1 double]
        X_g={};             %{q x 1} cells. Each cell c_i [n_i x 1 x T|1 x k double]
        Y_name={};          %{q x 1} cells of strings
        X_rg_name={};       %{q x 1} cells of strings
        X_beta_name={};     %{q x 1} cells. Each cell c_i {n_beta_i} cells of string
        X_time_name={};     %{q x 1} cells. Each cell c_i {n_time_i} cells of string
        X_g_name={};        %{q x k} cells. Each cell c_i {k} cells of strings
        simulated=0;        %1 if Y has been simulated
    end
    
    properties (SetAccess=private)
        N=[];                   %[1 x 1 integer] sum(n_1,...,n_q)
        T=[];                   %[1 x 1 integer] total number of time steps
        nvar=[];                %[1 x 1 integer] number of variables
        dim=[];                 %[q x 1 integer] number of rows of each Y
        %flag
        standardized=0;         %[1 x 1 boolean] 1 if Y has been standardized using the method standardize
        log_transformed=0;      %[1 x 1 boolean] 1 if Y has been log-transformed using the method log_transform
        X_rg_tv=1;              %[1 x 1 boolean] 1 if the remote loading are time variant 0 otherwise
        X_beta_tv=1;            %[1 x 1 boolean]
        X_time_tv=1;            %[1 x 1 boolean] 1 if the time loading are time variant 0 otherwise
        X_g_tv=1;               %[1 x 1 boolean] 1 if the ground loading are time variant 0 otherwise  
        %statistics
        Y_means={};             %{q x 1} cells. mean of the non-transfomed Y
        Y_stds={};              %{q x 1} cells. std of the non-transformed Y
        X_rg_means={};          %{q x 1} cells. mean of the non-transfomed X_rg
        X_rg_stds={};           %{q x 1} cells. std of the non-transformed X_rg
        X_beta_means={};        %{q x 1} cells. Each cell c_i [n_beta_i x 1 double] mean of the non-transfomed X_beta
        X_beta_stds={};         %{q x 1} cells. Each cell c_i [n_beta_i x 1 double] std of the non-transformed X_beta
        X_time_means={};        %{q x 1} cells. Each cell c_i [n_beta_i x 1 double] mean of the non-transfomed X_time
        X_time_stds={};         %{q x 1} cells. Each cell c_i [n_beta_i x 1 double] std of the non-transformed X_time
        X_g_means={};           %[q x k double] mean of the non-transfomed X_g
        X_g_stds={};            %[q x k double] std of the non-transformed X_g
    end
    
    methods
        function obj = stem_varset(Y,Y_name,X_rg,X_rg_name,X_beta,X_beta_name,X_time,X_time_name,X_g,X_g_name)
            if not(mod(nargin,2)==0)
                error('Not enough input arguments');
            end
            obj.Y=Y;
            obj.Y_name=Y_name;

            if nargin>=4
                if not(isempty(X_rg))
                    obj.X_rg=X_rg;
                    obj.X_rg_name=X_rg_name;
                end
            end
            if nargin>=6
                if not(isempty(X_beta))
                    obj.X_beta=X_beta;
                    obj.X_beta_name=X_beta_name;
                end
            end
            if nargin>=8
                if not(isempty(X_time))
                    obj.X_time=X_time;
                    obj.X_time_name=X_time_name;
                end
            end      
            if nargin>=10
                if not(isempty(X_g))
                    obj.X_g=X_g;
                    obj.X_g_name=X_g_name;
                end
            end              
        end

        function log_transform(obj)
            %only Y is log-transformed
            for i=1:length(obj.Y)
                temp=obj.Y{i};
                num1=sum(temp(:)<0);
                if num1>0
                    disp([num2str(num1),' negative value(s) are considered as zero']);
                    temp(temp(:)<0)=0;
                end
                num2=sum(temp(:)==0);
                if num2>0
                    disp([num2str(num2),' value(s) equal to zero are transformed to 0.05']);
                    temp(temp(:)==0)=0.05;
                end
                temp=log(temp);
                obj.Y{i}=temp;
            end
            obj.log_transformed=1;
        end
        
        function detrend(obj)
            for i=1:length(obj.Y)
                for j=1:size(obj.Y{i},1)
                    m1=nanmean(obj.Y{i}(j,:));
                    obj.Y{i}(j,:)=(obj.Y{i}(j,:)-m1);
                end
            end      
            for i=1:length(obj.X_beta)
                for j=1:size(obj.X_beta{i},2)
                    temp=squeeze(obj.X_beta{i}(:,j,:));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_beta{i}(:,j,:)=(obj.X_beta{i}(:,j,:)-m1)/std1;
                    obj.X_beta_means{i}(j)=m1;
                    obj.X_beta_stds{i}(j)=std1;
                end
            end           
            for i=1:length(obj.X_g)
                for j=1:size(obj.X_g{i},4)
                    temp=squeeze(obj.X_g{i}(:,:,:,j));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_g{i}(:,:,:,j)=(obj.X_g{i}(:,:,:,j)-m1)/std1;
                    obj.X_g_means{i}(j)=m1;
                    obj.X_g_stds{i}(j)=std1;
                end
            end                 
        end
        
        function standardize_sbs(obj)
            for i=1:length(obj.Y)
                m1=nanmean(obj.Y{i},2);
                std1=nanstd(obj.Y{i},1,2);
                for j=1:size(obj.Y{i},2)
                    obj.Y{i}(:,j)=(obj.Y{i}(:,j)-m1)./std1;
                end
            end
            
            for i=1:length(obj.X_beta)
                for j=1:size(obj.X_beta{i},2)
                    temp=squeeze(obj.X_beta{i}(:,j,:));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_beta{i}(:,j,:)=(obj.X_beta{i}(:,j,:)-m1)/std1;
                    obj.X_beta_means{i}(j)=m1;
                    obj.X_beta_stds{i}(j)=std1;
                end
            end          
            
            for i=1:length(obj.X_g)
                for j=1:size(obj.X_g{i},4)
                    temp=squeeze(obj.X_g{i}(:,:,:,j));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_g{i}(:,:,:,j)=(obj.X_g{i}(:,:,:,j)-m1)/std1;
                    obj.X_g_means{i}(j)=m1;
                    obj.X_g_stds{i}(j)=std1;
                end
            end                 
        end
        
        function standardize(obj)
            for i=1:length(obj.Y)
                m1=nanmean(obj.Y{i}(:));
                std1=nanstd(obj.Y{i}(:));
                obj.Y{i}=(obj.Y{i}-m1)/std1;
                obj.Y_means{i}=m1;
                obj.Y_stds{i}=std1;
            end
            for i=1:length(obj.X_rg)
                m1=mean(obj.X_rg{i}(:));
                std1=std(obj.X_rg{i}(:));
                if std1==0
                    m1=0;
                    std1=1;
                end
                obj.X_rg{i}=(obj.X_rg{i}-m1)/std1;
                obj.X_rg_means{i}=m1;
                obj.X_rg_stds{i}=std1;                
            end
            for i=1:length(obj.X_beta)
                for j=1:size(obj.X_beta{i},2)
                    temp=squeeze(obj.X_beta{i}(:,j,:));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_beta{i}(:,j,:)=(obj.X_beta{i}(:,j,:)-m1)/std1;
                    obj.X_beta_means{i}(j)=m1;
                    obj.X_beta_stds{i}(j)=std1;
                end
            end      
            
            for i=1:length(obj.X_time)
                for j=1:size(obj.X_time{i},2)
                    temp=squeeze(obj.X_time{i}(:,j,:));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_time{i}(:,j,:)=(obj.X_time{i}(:,j,:)-m1)/std1;
                    obj.X_time_means{i}(j)=m1;
                    obj.X_time_stds{i}(j)=std1;
                end
            end      
            
            for i=1:length(obj.X_g)
                for j=1:size(obj.X_g{i},4)
                    temp=squeeze(obj.X_g{i}(:,:,:,j));
                    m1=mean(temp(:));
                    std1=std(temp(:));
                    if std1==0
                        m1=0;
                        std1=1;
                    end
                    obj.X_g{i}(:,:,:,j)=(obj.X_g{i}(:,:,:,j)-m1)/std1;
                    obj.X_g_means{i}(j)=m1;
                    obj.X_g_stds{i}(j)=std1;
                end
            end               
            obj.standardized=1;
        end
        
        function index = get_Y_index(obj,name)
            index=find(strcmp(obj.Y_name,name));
        end

        %set methods
        
        function set.Y(obj,Y)
            if not(iscell(Y))
                error('Y must be a cell array');
            end
            for i=1:length(Y)
                if not(size(Y{i},2)==size(Y{1},2))
                    error('Each Y{i} must have the same number of temporal steps');
                end
                obj.dim(i)=size(Y{i},1);
            end
            obj.T=size(Y{1},2);
            obj.N=sum(obj.dim);
            obj.nvar=length(obj.dim);
            obj.Y=Y;
        end

        function set.Y_name(obj,Y_name)
            if not(iscell(Y_name))
                error('Y_name must be a cell array');
            end
            if not(length(Y_name)==length(obj.Y))
                error('The length of Y_name must be equal to length of Y');
            end
            obj.Y_name=Y_name;
        end
        
        function set.X_rg(obj,X_rg)
            if not(iscell(X_rg))
                error('X_rg must be a cell array');
            end
            if not(length(X_rg)==length(obj.Y))
                error('The number of cells of X_rg must be equal to the number of cells of Y');
            end
            for i=1:length(X_rg)
                if not(size(X_rg{i},1)==size(obj.Y{i},1))
                    error('X_rg{i} must have the same number of rows of Y{i}');
                end
                if not(size(X_rg{i},2)==1)
                    error('Each X_rg{i} must be a single covariate');
                end
                if not(size(X_rg{i},3)==obj.T || size(X_rg{i},3)==1)
                    error('Each X_rg{i} must have either 1 or T time steps');
                end
                if not(size(X_rg{1},3)==size(X_rg{i},3))
                    error('All the X_rg{i} must have the same temporal dimension');
                end
                if sum(isnan(X_rg{i}(:)))>0
                    error('X_rg cannot contain NaN');
                end
            end
            if size(X_rg{i},3)==1
                obj.X_rg_tv=0;
            end
            obj.X_rg=X_rg;
        end  
        
        function set.X_rg_name(obj,X_rg_name)
            if not(iscell(X_rg_name))
                error('X_rg_name must be a cell array');
            end
            if not(length(X_rg_name)==length(obj.X_rg))
                error('The length of X_rg_name must be equal to length of X_rg');
            end
            obj.X_rg_name=X_rg_name;
        end      
        
        function set.X_beta(obj,X_beta)
            if not(iscell(X_beta))
                error('X_beta must be a cell array');
            end
            if not(length(X_beta)==length(obj.Y))
                error('The number of cells of X_beta must be equal to the number of cells of Y');
            end
            for i=1:length(X_beta)
                if not(size(X_beta{i},1)==size(obj.Y{i},1))
                    error('X_beta{i} must have the same number of rows of Y{i}');
                end
                if not(size(X_beta{i},3)==obj.T || size(X_beta{i},3)==1)
                    error('Each X_beta{i} must have either 1 or T time steps');
                end
                if not(size(X_beta{1},3)==size(X_beta{i},3))
                    error('All the X_beta{i} must have the same temporal dimension');
                end
                if sum(isnan(X_beta{i}(:)))>0
                    error('X_beta cannot contain NaN');
                end                
            end
            if size(X_beta{1},3)==1
                obj.X_beta_tv=0;
            end
            obj.X_beta=X_beta;
        end  
        
        function set.X_beta_name(obj,X_beta_name)
            if not(iscell(X_beta_name))
                error('X_beta_name must be a cell array');
            end
            if not(length(X_beta_name)==length(obj.X_beta))
                error('The length of X_beta_name must be equal to length of X_beta');
            end
            for i=1:length(X_beta_name)
                if not(length(X_beta_name{i})==size(obj.X_beta{i},2))
                    error('The length of X_beta_name{i} must be equal to the number of covariates of X_beta{i}');
                end
            end
            obj.X_beta_name=X_beta_name;
        end      
        
        function set.X_time(obj,X_time)
            if not(iscell(X_time))
                error('X_time must be a cell array');
            end
            if not(length(X_time)==length(obj.Y))
                error('The number of cells of X_time must be equal to the number of cells of Y');
            end
            for i=1:length(X_time)
                if not(size(X_time{i},1)==size(obj.Y{i},1))
                    error('X_time{i} must have the same number of rows of Y{i}');
                end
                if not(size(X_time{i},3)==obj.T || size(X_time{i},3)==1)
                    error('Each X_time{i} must have either 1 or T time steps');
                end
                if not(size(X_time{i},3)==size(X_time{1},3))
                    error('All the X_time{i} must have the same temporal dimension');
                end
                if sum(isnan(X_time{i}(:)))>0
                    error('X_time cannot contain NaN');
                end                
            end
            if size(X_time{1},3)==1
                obj.X_time_tv=0;
            end
            obj.X_time=X_time;
        end  
        
        function set.X_time_name(obj,X_time_name)
            if not(iscell(X_time_name))
                error('X_time_name must be a cell array');
            end
            if not(length(X_time_name)==length(obj.X_time))
                error('The length of X_time_name must be equal to length of X_time');
            end
            for i=1:length(X_time_name)
                if not(length(X_time_name{i})==size(obj.X_time{i},2))
                    error('The length of X_time_name{i} must be equal to the number of covariates of X_time{i}');
                end
            end            
            obj.X_time_name=X_time_name;
        end      
        
        function set.X_g(obj,X_g)
            if not(iscell(X_g))
                error('X_g must be a cell array');
            end
            if not(length(X_g)==length(obj.Y))
                error('The number of cells of X_g must be equal to the number of cells of Y');
            end
            for i=1:length(X_g)
                if not(size(X_g{i},1)==size(obj.Y{i},1))
                    error('X_g{i} must have the same number of rows of Y{i}');
                end
                if not(size(X_g{i},2)==size(X_g{1},2))
                    error('Each X_g{i} must have the same number of covariates');
                end
                if not(size(X_g{i},3)==obj.T || size(X_g{i},3)==1)
                    error('Each X_g{i} must have either 1 or T time steps');
                end
                if not(size(X_g{i},3)==size(X_g{1},3))
                    error('All the X_g{i} must have the same temporal dimension');
                end
                if not(size(X_g{i},4)==size(X_g{1},4))
                    error('Each X_g{i} must have equal 4th dimension');
                end
                if sum(isnan(X_g{i}(:)))>0
                    error('X_g cannot contain NaN');
                end                
            end
            if size(X_g{1},3)==1
                obj.X_g_tv=0;
            end
            obj.X_g=X_g;
        end  
        
        function set.X_g_name(obj,X_g_name)
            if not(iscell(X_g_name))
                error('X_g_name must be a cell array');
            end
            if not(length(X_g_name)==length(obj.X_g))
                error('The length of X_g_name must be equal to length of X_g');
            end
            for i=1:length(X_g_name)
                if not(length(X_g_name{i})==size(obj.X_g{i},4))
                    error('The length of X_g_name{i} must be equal to k');
                end
            end
            obj.X_g_name=X_g_name;
        end           
        
    end
end