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
    
    %CONSTANTS
    %ni   - the number of sites for the i-th variable, i=1,...,q
    %ni_b - the number of loading vectors for the i-th variable related to the beta parameter
    %ni_t - the number of loading vectors for the i-th variable related to the latent variable z
    %K    - the number of loading vectors related to the latent variable w_g
    %N    - n1+...+nq total number of observation sites for all the variables
    %T - number of temporal steps
    %TT = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    %p - dimension of the latent temporal variable z    
    
    properties
        Y={};               %[double]   {qx1}(nixT)         observed data 
        X_rg={};            %[double]   {qx1}(nix1xTT)      loading vectors related to the latent variable w_r 
        X_beta={};          %[double]   {qx1}(nixni_bxTT)   loading vectors related to the beta parameter
        X_time={};          %[double]   {qx1}(nixni_txTT)   loading vectors related to the latent variable z
        X_g={};             %[double]   {qx1}(nix1xTTxK)    loading vectors related to the latent variable w_g
        Y_name={};          %[string]   {qx1}               variable names
        X_rg_name={};       %[string]   {qx1}               name of the loading vectors related to the latent variable w_r
        X_beta_name={};     %[string]   {qxni_b}            name of the loading vectors related to the beta parameter
        X_time_name={};     %[string]   {qxni_t}            name of the loading vectors related to the latent variable z
        X_g_name={};        %[string]   {qxK}               name of the loading vectors related to the latent variable w_g
        simulated=0;        %[boolean]  (1x1)               1: the Y data are simulated; 0: the Y data are observed data
    end
    
    properties (SetAccess=private)
        N=[];                   %[integer >0]      (1x1) n1+...+nq
        T=[];                   %[integer >0]      (1x1) total number of time steps
        nvar=[];                %[integer >0]      (1x1) total number of variables (q)
        dim=[];                 %[integer]         (qx1) number of time series for each variable

        standardized=0;         %[boolean]         (1x1) 1: Y, X_rg, X_beta, X_time and X_g has been standardized; 0: otherwise
        log_transformed=0;      %[boolean]         (1x1) 1: Y has been log-transformed using the method log_transform; 0: otherwise
        X_rg_tv=1;              %[boolean]         (1x1) 1: the loading vectors related to the latent variable w_r are time variant; 0:otherwise
        X_beta_tv=1;            %[boolean]         (1x1) 1: the loading vectors related to the beta parameter are time variant; 0:otherwise
        X_time_tv=1;            %[boolean]         (1x1) 1: the loading vectors related to the latent variable z are time variant; 0:otherwise
        X_g_tv=1;               %[boolean]         (1x1) 1: the loading vectors related to the latent variable w_g are time variant; 0:otherwise
        
        Y_means={};             %[double]          {qx1} averages of the non-standardized Y
        Y_stds={};              %[double]          {qx1} standard deviations of the non-standardized Y
        X_rg_means={};          %[double]          {qx1} averages of the non-standardized X_rg
        X_rg_stds={};           %[double]          {qx1} standard deviations of the non-standardized X_rg
        X_beta_means={};        %[double]          {qx1}(ni_bx1) averages of the non-standardized X_beta
        X_beta_stds={};         %[double]          {qx1}(ni_bx1) standard deviations of the non-standardized X_beta
        X_time_means={};        %[double]          {qx1}(ni_tx1) averages of the non-standardized X_time
        X_time_stds={};         %[double]          {qx1}(ni_tx1) standard deviations of the non-standardized X_time
        X_g_means={};           %[double]          (qxK) averages of the non-standardized X_g
        X_g_stds={};            %[double]          {qx1}(ni_tx1) standard deviations of the non-standardized X_g
    end
    
    methods
        function obj = stem_varset(Y,Y_name,X_rg,X_rg_name,X_beta,X_beta_name,X_time,X_time_name,X_g,X_g_name)
            %DESCRIPTION: is the constructor of the class stem_varset
            %
            %INPUT
            %            
            %Y               -  [double]   {qx1}(nixT)         observed data
            %Y_name          -  [string]   {qx1}               variable names
            %X_rg            -  [double]   {qx1}(nix1xTT)      loading vectors related to the latent variable w_r
            %X_rg_name       -  [string]   {qx1}               name of the loading vectors related to the latent variable w_r
            %X_beta          -  [double]   {qx1}(nixni_bxTT)   loading vectors related to the beta parameter
            %X_beta_name     -  [string]   {qxni_b}            name of the loading vectors related to the beta parameter
            %X_time          -  [double]   {qx1}(nixni_txTT)   loading vectors related to the latent variable z
            %X_time_name     -  [string]   {qxni_t}            name of the loading vectors related to the latent variable z
            %X_g             -  [double]   {qx1}(nix1xTTxK)    loading vectors related to the latent variable w_g
            %X_g_name        -  [string]   {qxK}               name of the loading vectors related to the latent variable w_g
            %
            %OUTPUT
            %obj             - [stem_varset object] (1x1) the stem_varset object
            
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
            %DESCRIPTION: log-transforms the matrix Y
            %
            %INPUT
            %obj - [stem_varset object] (1x1) the stem_varset object
            %
            %OUTPUT
            %
            %none: the matrix Y is updated              

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
        
        function detrend_Y(obj)
            %DESCRIPTION: remove the mean from each time series in Y
            %
            %INPUT
            %obj - [stem_varset object] (1x1) the stem_varset object
            %
            %OUTPUT
            %
            %none: the Y property is updated            
            for i=1:length(obj.Y)
                for j=1:size(obj.Y{i},1)
                    m1=nanmean(obj.Y{i}(j,:));
                    obj.Y{i}(j,:)=(obj.Y{i}(j,:)-m1);
                end
            end                    
        end
        
        function standardize_Y(obj)
            %DESCRIPTION: each time series in Y is standardized
            %
            %INPUT
            %obj - [stem_varset object] (1x1) the stem_varset object
            %
            %OUTPUT
            %
            %none: the Y property is updated        
            
            for i=1:length(obj.Y)
                m1=nanmean(obj.Y{i},2);
                std1=nanstd(obj.Y{i},1,2);
                for j=1:size(obj.Y{i},2)
                    obj.Y{i}(:,j)=(obj.Y{i}(:,j)-m1)./std1;
                end
            end                
        end
        
        function standardize(obj)
            %DESCRIPTION: standardize the matrices Y, X_rg, X_beta, X_time and X_g with respect to their overall mean and overall standard deviation
            %
            %INPUT
            %obj - [stem_varset object] (1x1) the stem_varset object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
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
            %DESCRIPTION: returns the index of the variable given its name
            %
            %INPUT
            %obj    - [stem_varset object]  (1x1) the stem_varset object
            %name   - [string]              (1x1) the name of the variable
            %
            %OUTPUT
            %
            %none: the Y property is updated     
            
            index=find(strcmp(obj.Y_name,name));
        end

        %Class set methods
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