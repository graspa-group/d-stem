%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
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
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_varset < handle
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    %
    %CONSTANTS
    %q     - the number of variables
    %n_i   - the number of sites for the i-th variable, i=1,...,q
    %nb_i  - the number of loading vectors for the i-th variable related to the beta parameter
    %nz_i  - the number of loading vectors for the i-th variable related to the latent variable z
    %nw    - the number of latent variables w
    %N     - N=n_1+...+n_q total number of observation sites for all the variables
    %T     - number of time steps
    %TT    - TT=T if the loading coefficients are time-variant and TT=1 if they are time-invariant

    properties
        Y={};               %[double]   {q}(n_i x T)           observed data 
        X_bp={};            %[double]   {q}(n_i x TT)          loading vectors related to the latent variable w_b 
        X_beta={};          %[double]   {q}(n_i x nb_i x TT)   loading vectors related to the beta parameter
        X_z={};             %[double]   {q}(n_i x nz_i x TT)   loading vectors related to the latent variable z
        X_p={};             %[double]   {q}(n_i x nw x TT)     loading vectors related to the latent variable w_p
        X_h={};             %[double]   {q}(n_i x TT)          domain of the functional observations in Y 
        Y_name={};          %[string]   {q}                    variable names
        X_bp_name={};       %[string]   {q}                    name of the loading vectors related to the latent variable w_b
        X_beta_name={};     %[string]   {q}{nb_i}              name of the loading vectors related to the beta parameter
        X_z_name={};        %[string]   {q}{nz_i}              name of the loading vectors related to the latent variable z
        X_p_name={};        %[string]   {q}{nw}                name of the loading vectors related to the latent variable w_p        
        X_h_name=[];        %[string]   (1x1)                  name of the domain of the function observations in Y
        
        Y_unit={};          %[string]   (1x1)                  unit of variable y
        X_beta_unit={};     %[string]   (1x1)                  unit of covariates
        X_h_unit=[];        %[string]   (1x1)                  unit of function domain 
        T_unit=[];          %[string]   (1x1)                  unit of time step

        simulated=0;        %[boolean]  (1x1)                  1: the Y data are simulated; 0: the Y data are observed data
    end
    
    properties (SetAccess=private)
        standardized=0;         %[boolean]         (1x1) 1: Y, X_bp, X_beta, X_z and X_p has been standardized; 0: otherwise
        log_transformed=0;      %[boolean]         (1x1) 1: Y has been log-transformed using the method log_transform; 0: otherwise
        boxcox_transformed=0;   %[boolean]         (1x1) 1: Y has been boxcox-transformed using the method boxcox_transform; 0: otherwise
        
        Y_means={};             %[double]          {q} averages of the non-standardized Y
        Y_stds={};              %[double]          {q} standard deviations of the non-standardized Y
        Y_lambda={};            %[double]          {q} the lambda parameters of the boxcox transformation
        X_bp_means={};          %[double]          {q} averages of the non-standardized X_bp
        X_bp_stds={};           %[double]          {q} standard deviations of the non-standardized X_bp
        X_beta_means={};        %[double]          {q}(1 x nb_i) averages of the non-standardized X_beta
        X_beta_stds={};         %[double]          {q}(1 x nb_i) standard deviations of the non-standardized X_beta
        X_z_means={};           %[double]          {q}(1 x nz_i) averages of the non-standardized X_z
        X_z_stds={};            %[double]          {q}(1 x nz_i) standard deviations of the non-standardized X_z
        X_p_means={};           %[double]          {q}(1 x nw) averages of the non-standardized X_p
        X_p_stds={};            %[double]          {q}(1 x nw) standard deviations of the non-standardized X_p
    end
    
    properties (Dependent, SetAccess = private)
        nvar=[];                %[integer >0]      (1x1) total number of variables (q)
        dim=[];                 %[integer]         (qx1) number of spatial sites for each variable
        N=[];                   %[integer >0]      (1x1) n_1+...+n_q
        T=[];                   %[integer >0]      (1x1) total number of time steps
        
        X_bp_tv=[];             %[boolean]         (1x1) 1: the loading vectors related to the latent variable w_b are time variant; 0:otherwise
        X_beta_tv=[];           %[boolean]         (1x1) 1: the loading vectors related to the beta parameter are time variant; 0:otherwise
        X_z_tv=[];              %[boolean]         (1x1) 1: the loading vectors related to the latent variable z are time variant; 0:otherwise
        X_p_tv=[];              %[boolean]         (1x1) 1: the loading vectors related to the latent variable w_p are time variant; 0:otherwise
        X_h_tv=[];              %[boolean]         (1x1) 1: the domain of the functional observartions is time variant; 0:otherwise
    end
    
    methods
        function obj = stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name,X_h,X_h_name)
            %DESCRIPTION: is the constructor of the class stem_varset
            %
            %INPUT
            %
            %Y               -  [double]   {q}(n_i x T)          observed data
            %Y_name          -  [string]   {q}                   variable names
            %X_bp            -  [double]   {q}(n_i x TT)         loading vectors related to the latent variable w_b
            %X_bp_name       -  [string]   {q}                   name of the loading vectors related to the latent variable w_b
            %X_beta          -  [double]   {q}(n_i x nb_i x TT)  loading vectors related to the beta parameter
            %X_beta_name     -  [string]   {q}{nb_i}             name of the loading vectors related to the beta parameter
            %X_z             -  [double]   {q}(n_i x nz_i x TT)  loading vectors related to the latent variable z
            %X_z_name        -  [string]   {q}{nz_i}             name of the loading vectors related to the latent variable z
            %X_p             -  [double]   {q}(n_i x nw x TT)    loading vectors related to the latent variable w_p
            %X_p_name        -  [string]   {q}{nw}               name of the loading vectors related to the latent variable w_p
            %X_h             -  [double]   {q}(n_i x TT)         domain of the functional observations in Y 
            %X_h_name        -  [string]   1x1                   name of the domain of the function observations in Y
            %
            %OUTPUT
            %obj             -  [stem_varset object] (1x1) the stem_varset object
            
            temp_struct=Y;
            if isstruct(temp_struct)
                names=fieldnames(temp_struct);
                if sum(strcmp('Y',names))
                    Y=temp_struct.Y;
                end
                if sum(strcmp('Y_name',names))
                    Y_name=temp_struct.Y_name;
                end
                if sum(strcmp('X_bp',names))
                    X_bp=temp_struct.X_bp;
                else
                    X_bp=[];
                end 
                if sum(strcmp('X_bp_name',names))
                    X_bp_name=temp_struct.X_bp_name;
                else
                    X_bp_name=[];
                end 
                if sum(strcmp('X_beta',names))
                    X_beta=temp_struct.X_beta;
                else
                    X_beta=[];
                end 
                if sum(strcmp('X_beta_name',names))
                    X_beta_name=temp_struct.X_beta_name;
                else
                    X_beta_name=[];
                end 
                if sum(strcmp('X_z',names))
                    X_z=temp_struct.X_z;
                else
                    X_z=[];
                end 
                if sum(strcmp('X_z_name',names))
                    X_z_name=temp_struct.X_z_name;
                else
                    X_z_name=[];
                end 
                if sum(strcmp('X_p',names))
                    X_p=temp_struct.X_p;
                else
                    X_p=[];
                end 
                if sum(strcmp('X_p_name',names))
                    X_p_name=temp_struct.X_p_name;
                else
                    X_p_name=[];
                end 
                if sum(strcmp('X_h',names))
                    X_h=temp_struct.X_h;
                else
                    X_h=[];
                end 
                if sum(strcmp('X_h_name',names))
                    X_h_name=temp_struct.X_h_name;
                else
                    X_h_name=[];
                end 
            else
                if not(mod(nargin,2)==0)
                    error('Not enough input arguments');
                end
            end
            
            for i=1:length(Y)
                if not(size(Y{i},2)==size(Y{1},2))
                    error('Each Y{i} must have the same number of temporal steps');
                end
            end
            obj.Y=Y;
            
            if not(length(Y_name)==length(Y))
                error('The length of Y_name must be equal to length of Y');
            end
            obj.Y_name=Y_name;
            
            if nargin>2||isstruct(temp_struct)
                if not(isempty(X_bp))
                    if not(length(X_bp)==length(Y))
                        error('The number of cells of X_bp must be equal to the number of cells of Y');
                    end
                    for i=1:length(X_bp)
                        if not(size(X_bp{i},1)==size(Y{i},1))
                            error('X_bp{i} must have the same number of rows of Y{i}');
                        end
                        if not(size(X_bp{i},2)==obj.T || size(X_bp{i},2)==1)
                            error('Each X_bp{i} must have either 1 or T time steps');
                        end
                        if not(size(X_bp{1},2)==size(X_bp{i},2))
                            error('All the X_bp{i} must have the same number of columns');
                        end
                        if sum(isnan(X_bp{i}(:)))>0
                            error('X_bp cannot contain NaN');
                        end
                    end
                    obj.X_bp=X_bp;
                    
                    if not(length(X_bp_name)==length(X_bp))
                        error('The length of X_bp_name must be equal to length of X_bp');
                    end
                    obj.X_bp_name=X_bp_name;
                end
            end
            
            if nargin>4||isstruct(temp_struct)
                if not(isempty(X_beta))
                    if not(length(X_beta)==length(Y))
                        error('The number of cells of X_beta must be equal to the number of cells of Y');
                    end
                    for i=1:length(X_beta)
                        T_max=[];
                        if not(isempty(X_beta{i}))
                            T_max=size(X_beta{i},3);
                        end
                    end

                    for i=1:length(X_beta)
                        if not(isempty(X_beta{i}))
                            if not(size(X_beta{i},1)==size(Y{i},1))
                                error('X_beta{i} must have the same number of rows of Y{i}');
                            end
                            if not(size(X_beta{i},3)==obj.T || size(X_beta{i},3)==1)
                                error('Each X_beta{i} must have either 1 or T time steps');
                            end
                            if not(size(X_beta{i},3)==T_max)
                                error('All the X_beta{i} must have the same temporal dimension');
                            end
                        else
                            if not(isempty(T_max))
                                X_beta{i}=zeros(size(Y{i},1),1,T_max);
                            end
                        end
                    end
                    obj.X_beta=X_beta;
                    
                    if not(length(X_beta_name)==length(X_beta))
                        error('The length of X_beta_name must be equal to length of X_beta');
                    end
                    for i=1:length(X_beta_name)
                        if not(isempty(X_beta_name{i}))
                            if not(length(X_beta_name{i})==size(X_beta{i},2))
                                error('The length of X_beta_name{i} must be equal to the number of covariates of X_beta{i}');
                            end
                        end
                    end
                    obj.X_beta_name=X_beta_name;
                end
            end
            
            if nargin>6||isstruct(temp_struct)
                if not(isempty(X_z))
                    if not(length(X_z)==length(Y))
                        error('The number of cells of X_z must be equal to the number of cells of Y');
                    end
                    for i=1:length(X_z)
                        T_max=[];
                        if not(isempty(X_z{i}))
                            T_max=size(X_z{i},3);
                        end
                    end
                    for i=1:length(X_z)
                        if not(isempty(X_z{i}))
                            if not(size(X_z{i},1)==size(Y{i},1))
                               error('X_z{i} must have the same number of rows of Y{i}');
                            end
                            if not(size(X_z{i},3)==obj.T || size(X_z{i},3)==1)
                                error('Each X_z{i} must have either 1 or T time steps');
                            end
                            if not(size(X_z{i},3)==size(X_z{1},3))
                                error('All the X_z{i} must have the same temporal dimension');
                            end
                            if sum(isnan(X_z{i}(:)))>0
                                error('X_z cannot contain NaN');
                            end
                        else
                            X_z{i}=zeros(size(Y{i},1),1,T_max);
                        end
                    end
                    obj.X_z=X_z;
                    
                    if not(length(X_z_name)==length(X_z))
                        error('The length of X_z_name must be equal to length of X_z');
                    end
                    for i=1:length(X_z_name)
                        if not(isempty(X_z_name{i}))
                            if not(length(X_z_name{i})==size(X_z{i},2))
                                error('The length of X_z_name{i} must be equal to the number of covariates of X_z{i}');
                            end
                        end
                    end
                    obj.X_z_name=X_z_name;
                end
            end
            
            if nargin>8||isstruct(temp_struct)
                if not(isempty(X_p))
                    if not(length(X_p)==length(Y))
                        error('The number of cells of X_p must be equal to the number of cells of Y');
                    end
                    for i=1:length(X_p)
                        if not(size(X_p{i},1)==size(Y{i},1))
                            error('X_p{i} must have the same number of rows of Y{i}');
                        end
                        if not(size(X_p{i},2)==size(X_p{1},2))
                            error('Each X_p{i} must have the same number of columns');
                        end
                        if not(size(X_p{i},3)==obj.T || size(X_p{i},3)==1)
                            error('Each X_p{i} must have either 1 or T time steps');
                        end
                        if not(size(X_p{i},3)==size(X_p{1},3))
                            error('All the X_p{i} must have the same temporal dimension');
                        end
                        if sum(isnan(X_p{i}(:)))>0
                            error('X_p cannot contain NaN');
                        end
                    end
                    obj.X_p=X_p;
                    
                    if not(length(X_p_name)==length(X_p))
                        error('The length of X_p_name must be equal to length of X_p');
                    end
                    for i=1:length(X_p_name)
                        if not(length(X_p_name{i})==size(X_p{i},2))
                            error(['The length of X_p_name{i} must be equal to ',num2str(size(X_p{i},2))]);
                        end
                    end
                    obj.X_p_name=X_p_name;
                end
            end
            
            if nargin>10||isstruct(temp_struct)
                if not(isempty(X_h))
                    if not(length(X_h)==length(Y))
                        error('The number of cells of X_h must be equal to the number of cells of Y');
                    end
                    for i=1:length(X_h)
                        if not(size(X_h{i},1)==size(Y{i},1))
                            error('X_h{i} must have the same number of rows of Y{i}');
                        end
                        if not(size(X_h{i},2)==obj.T || size(X_h{i},2)==1)
                            error('Each X_h{i} must have either 1 or T time steps');
                        end
                    end
                    obj.X_h=X_h;
                    
                    if not(ischar(X_h_name))
                        error('X_h_name must be a string');
                    end
                    obj.X_h_name=X_h_name;
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
                num1=sum(temp(:)<=0);
                if num1>0
                    disp([num2str(num1),' null and negative value(s) are considered as NaN']);
                    temp(temp(:)<=0)=NaN;
                end
                temp=log(temp);
                obj.Y{i}=temp;
            end
            obj.log_transformed=1;
        end
        
        function boxcox_transform(obj,lambda)
            %DESCRIPTION: boxcox-transforms the matrix Y
            %
            %INPUT
            %obj    - [stem_varset object] (1x1) the stem_varset object
            %lambda - [double]             (qx1) the vector of lambda paramters
            %
            %OUTPUT
            %
            %none: the matrix Y is updated              

            for i=1:length(obj.Y)
                if isnan(lambda(i))
                    temp=obj.Y{i};
                    temp(temp<=0)=NaN;
                    temp(isnan(temp))=[];
                    [~,lambda_estimated]=boxcox(temp(:));
                else
                    lambda_estimated=lambda(i);
                end
                obj.Y_lambda{i}=lambda_estimated;
                
                temp=obj.Y{i};
                temp(temp<=0)=NaN;
                temp=(temp.^lambda_estimated-1)./lambda_estimated;
               
                obj.Y{i}=temp;
            end
            obj.boxcox_transformed=1;
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
                obj.Y_means{i}=m1;
                obj.Y_stds{i}=std1;
            end    
        end
        
        function standardize(obj)
            %DESCRIPTION: standardize the matrices Y, X_bp, X_beta, X_z and X_p with respect to their overall mean and overall standard deviation
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
            for i=1:length(obj.X_bp)
                m1=mean(obj.X_bp{i}(:));
                std1=std(obj.X_bp{i}(:));
                if std1==0
                    m1=0;
                    std1=1;
                end
                obj.X_bp{i}=(obj.X_bp{i}-m1)/std1;
                obj.X_bp_means{i}=m1;
                obj.X_bp_stds{i}=std1;                
            end
            for i=1:length(obj.X_beta)
                for j=1:size(obj.X_beta{i},2)
                    temp=squeeze(obj.X_beta{i}(:,j,:));
                    if sum(not(temp==1)&not(temp==0))==0
                        %if the matrix includes only 0 and 1 do not standardize
                        m1=0;
                        std1=1;
                    else
                        m1=nanmean(temp(:));
                        std1=nanstd(temp(:));
                        if std1==0
                            m1=0;
                            std1=1;
                        end
                    end
                    obj.X_beta{i}(:,j,:)=(obj.X_beta{i}(:,j,:)-m1)/std1;
                    obj.X_beta_means{i}(j)=m1;
                    obj.X_beta_stds{i}(j)=std1;
                end
            end      
            
            for i=1:length(obj.X_z)
                for j=1:size(obj.X_z{i},2)
                    temp=squeeze(obj.X_z{i}(:,j,:));
                    if sum(not(temp==1)&not(temp==0))==0
                        %if the matrix includes only 0 and 1 do not standardize
                        m1=0;
                        std1=1;
                    else
                        m1=mean(temp(:));
                        std1=std(temp(:));
                        if std1==0
                            m1=0;
                            std1=1;
                        end
                    end
                    obj.X_z{i}(:,j,:)=(obj.X_z{i}(:,j,:)-m1)/std1;
                    obj.X_z_means{i}(j)=m1;
                    obj.X_z_stds{i}(j)=std1;
                end
            end      
            
            for i=1:length(obj.X_p)
                for j=1:size(obj.X_p{i},2)
                    temp=squeeze(obj.X_p{i}(:,j,:));
                    if sum(not(temp==1)&not(temp==0))==0
                        %if the matrix includes only 0 and 1 do not standardize
                        m1=0;
                        std1=1;
                    else
                        m1=mean(temp(:));
                        std1=std(temp(:));
                        if std1==0
                            m1=0;
                            std1=1;
                        end
                    end
                    obj.X_p{i}(:,j,:)=(obj.X_p{i}(:,j,:)-m1)/std1;
                    obj.X_p_means{i}(j)=m1;
                    obj.X_p_stds{i}(j)=std1;
                end
            end               
            obj.standardized=1;
        end
        
        function time_crop(obj,indices)
            %DESCRIPTION: crop the matrices Y, X_bp, X_beta, X_z, X_p and X_h with respect to time
            %
            %INPUT
            %obj         - [stem_varset object]       (1x1) the stem_data object
            %indices     - [integer >0]               (dTx1) a dTx1 vector of (possibly non consecutive) temporal indices
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated            
            for i=1:length(obj.Y)
                obj.Y{i}=obj.Y{i}(:,indices);
            end
            if not(isempty(obj.X_beta))
                if obj.X_beta_tv
                    for i=1:length(obj.X_beta)
                        obj.X_beta{i}=obj.X_beta{i}(:,:,indices);
                    end
                end
            end
            if not(isempty(obj.X_bp))
                if obj.X_bp_tv
                    for i=1:length(obj.X_bp)
                        obj.X_bp{i}=obj.X_bp{i}(:,:,indices);
                    end
                end
            end               
            if not(isempty(obj.X_z))
                if obj.X_z_tv
                    for i=1:length(obj.X_z)
                        obj.X_z{i}=obj.X_z{i}(:,:,indices);
                    end
                end
            end
            if not(isempty(obj.X_p))
                if obj.X_p_tv
                    for i=1:length(obj.X_p)
                        obj.X_p{i}=obj.X_p{i}(:,:,indices);
                    end
                end
            end
            if not(isempty(obj.X_h))
                if obj.X_h_tv
                    for i=1:length(obj.X_h)
                        obj.X_h{i}=obj.X_h{i}(:,indices);
                    end
                end
            end
        end     
        
        function time_average(obj,n_steps)
            %DESCRIPTION: computes time averages of n_steps for the matrice the matrices Y, X_bp, X_beta, X_z, X_p and X_h
            %
            %INPUT
            %obj        - [stem_varset object] (1x1) the stem_varset object
            %n_steps    - [integer >0]         (1x1) the number of temporal steps to average
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            indices=0:n_steps:obj.T;
            if not(indices(end)==obj.T)
                disp(['The last ',num2str(obj.T-indices(end)),' time steps will be discarded when taking time averages']);
            end

            disp('Time averaging started...');
            Y_temp=cell(length(obj.Y),1);
            for i=1:length(obj.Y)
                for j=1:length(indices)-1
                    Y_temp{i}(:,j)=nanmean(obj.Y{i}(:,indices(j)+1:indices(j+1)),2);
                end
            end
            obj.Y=Y_temp;

            if not(isempty(obj.X_bp))
                if obj.X_bp_tv
                    X_bp_temp=cell(length(obj.X_bp),1);
                    for i=1:length(obj.X_bp)
                        for j=1:length(indices)-1
                            X_bp_temp{i}(:,:,j)=nanmean(obj.X_bp{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.X_bp=X_bp_temp;
                end
            end
            
            if not(isempty(obj.X_beta))
                if obj.X_beta_tv
                    X_beta_temp=cell(length(obj.X_beta),1);
                    for i=1:length(obj.X_beta)
                        for j=1:length(indices)-1
                            X_beta_temp{i}(:,:,j)=nanmean(obj.X_beta{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.X_beta=X_beta_temp;
                end
            end
            
            if not(isempty(obj.X_z))
                if obj.X_z_tv
                    X_z_temp=cell(length(obj.X_z),1);
                    for i=1:length(obj.X_z)
                        for j=1:length(indices)-1
                            X_z_temp{i}(:,:,j)=nanmean(obj.X_z{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.X_z=X_z_temp;
                end
            end
            
            if not(isempty(obj.X_p))
                if obj.X_p_tv
                    X_p_temp=cell(length(obj.X_p),1);
                    for i=1:length(obj.X_p)
                        for j=1:length(indices)-1
                            X_p_temp{i}(:,:,j)=nanmean(obj.X_p{i}(:,:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.X_p=X_p_temp;
                end
            end
            
            if not(isempty(obj.X_h))
                if obj.X_h_tv
                    X_h_temp=cell(length(obj.X_h),1);
                    for i=1:length(obj.X_h)
                        for j=1:length(indices)-1
                            X_h_temp{i}(:,j)=nanmean(obj.X_h{i}(:,indices(j)+1:indices(j+1)),3);
                        end
                    end
                    obj.X_h=X_h_temp;
                end
            end
            disp('Time averaging ended.');
        end        
         
        function site_crop(obj,var_name,indices)
            %DESCRIPTION: remove specific sites from the dataset
            %
            %INPUT
            %obj                 - [stem_varset object] (1x1) the stem_varset object
            %var_name            - [string]             (1x1) the name of the variable from which to remove the sites
            %indices             - [integer >0]         (dNx1) the indices of the sites to remove
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z, X_p and X_h are updated
            
            if not(isempty(obj.X_h_name)) %f-HDGM
                idx_var=obj.get_Y_index(var_name)';
                if isempty(idx_var)
                    error('Variable not found');
                end
                if max(indices>obj.N)
                    error(['The maximum value of indices cannot be greater than ',num2str(obj.N)]);
                end
                for idx=idx_var
                    obj.Y{idx}(indices,:)=[];
                    if not(isempty(obj.X_bp))
                        obj.X_bp{idx}(indices,:,:)=[];
                    end
                    if not(isempty(obj.X_beta))
                        obj.X_beta{idx}(indices,:,:)=[];
                    end
                    if not(isempty(obj.X_z))
                        obj.X_z{idx}(indices,:,:)=[];
                    end
                    if not(isempty(obj.X_p))
                        obj.X_p{idx}(indices,:,:)=[];
                    end
                    if not(isempty(obj.X_h))
                        obj.X_h{idx}(indices,:)=[];
                    end
                end
            else
                idx_var=obj.get_Y_index(var_name);
                if isempty(idx_var)
                    error('Variable not found');
                end
                if max(indices>obj.N)
                    error(['The maximum value of indices cannot be greater than ',num2str(obj.N)]);
                end
                obj.Y{idx_var}(indices,:)=[];

                if not(isempty(obj.X_bp))
                    obj.X_bp{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.X_beta))
                    obj.X_beta{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.X_z))
                    obj.X_z{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.X_p))
                    obj.X_p{idx_var}(indices,:,:)=[];
                end
                if not(isempty(obj.X_h))
                    obj.X_h{idx_var}(indices,:)=[];
                end
            end
            
        end
        
        function idx_vec = missing_crop(obj,threshold)
            %DESCRIPTION: remove sites with a missing data rate higher or equal to a threshold
            %
            %INPUT
            %obj                 - [stem_varset object]   (1x1) the stem_varset object
            %threshold           - [double >0 and <1]     (1x1) the threshold
            %
            %OUTPUT
            %
            %idx_vec             - [integer >0]           {qx1} the cell-array with the indices of the removed sites 
            idx_vec=cell(obj.nvar,1);
            for i=1:obj.nvar
                idx_vec{i}=find((sum(isnan(obj.Y{i}),2)/obj.T)>=threshold);
                if not(isempty(idx_vec{i}))
                    obj.site_crop(obj.Y_name{i},idx_vec{i});
                    disp(['    Removed ',num2str(length(idx_vec{i})),' of ',num2str(obj.dim(i)),' sites for the variable ',obj.Y_name{i}]);
                end
            end
        end
           
        function permute(obj,indices)
            %DESCRIPTION: permute rows of Y, X_bp, X_beta, X_z, X_p and X_h with respect to indices
            %
            %INPUT
            %obj         - [stem_varset object]       (1x1) the stem_data object
            %indices     - [integer >0]               {q}(n_ix1)|(n_ix1) a cell array of indices or o vector of indices. In the latter case, the same vector is applied to all the q variables
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated        
            if iscell(indices)
                if not(length(indices))==obj.nvar
                    error(['The cell array indices must have ',num2str(obj.nvar),' cells']);
                end
                dims=obj.dim;
                for i=1:length(indices)
                    if not(length(indices{i})==dims(i))
                        error('The length of each vector in the cell array indices must be equal to the number of observations for each variable');
                    end
                    if min(indices{i})<1
                        error('The elements of each vector in the cell array indices must be > 0');
                    end
                    if max(indices{i})>dims(i)
                        error('The elements of each vector in the cell array indices cannot exceed the number of observations for each variable');
                    end
                end
            else
                dims=obj.dim;
                if min(indices)<1
                    error('The elements of indices must be > 0');
                end
                for i=1:length(dims)
                    if not(length(indices)==dims(i))
                        error(['The length of indices is not comformable to the number of observations in variable number ',num2str(i)]);
                    end
                    if max(indices)>dims(i)
                        error('The elements of indices cannot exceed the number of observations for each variable');
                    end
                end
            end
            
            for i=1:length(obj.Y)
                if iscell(indices)
                    obj.Y{i}=obj.Y{i}(indices{i},:);
                else
                    obj.Y{i}=obj.Y{i}(indices,:);
                end
            end
            if not(isempty(obj.X_beta))
                if obj.X_beta_tv
                    for i=1:length(obj.X_beta)
                        if iscell(indices)
                            obj.X_beta{i}=obj.X_beta{i}(indices{i},:,:);
                        else
                            obj.X_beta{i}=obj.X_beta{i}(indices,:,:);
                        end
                    end
                end
            end
            if not(isempty(obj.X_bp))
                if obj.X_bp_tv
                    for i=1:length(obj.X_bp)
                        if iscell(indices)
                            obj.X_bp{i}=obj.X_bp{i}(indices{i},:,:);
                        else
                            obj.X_bp{i}=obj.X_bp{i}(indices,:,:);
                        end
                    end
                end
            end               
            if not(isempty(obj.X_z))
                if obj.X_z_tv
                    for i=1:length(obj.X_z)
                        if iscell(indices)
                            obj.X_z{i}=obj.X_z{i}(indices{i},:,:);
                        else
                            obj.X_z{i}=obj.X_z{i}(indices,:,:);
                        end
                    end
                end
            end
            if not(isempty(obj.X_p))
                if obj.X_p_tv
                    for i=1:length(obj.X_p)
                        if iscell(indices)
                            obj.X_p{i}=obj.X_p{i}(indices{i},:,:);
                        else
                            obj.X_p{i}=obj.X_p{i}(indices,:,:);
                        end
                    end
                end
            end
            if not(isempty(obj.X_h))
                if obj.X_h_tv
                    for i=1:length(obj.X_h)
                        if iscell(indices)
                            obj.X_h{i}=obj.X_h{i}(indices{i},:);
                        else
                            obj.X_h{i}=obj.X_h{i}(indices,:);
                        end
                    end
                end
            end
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
            %none: the index of the requested variable     
            
            index=find(strcmp(obj.Y_name,name));
        end
        
        function index = get_X_beta_index(obj,name,variable_idx)
            %DESCRIPTION: returns the index of the variable given its name
            %
            %INPUT
            %obj          - [stem_varset object]  (1x1) the stem_varset object
            %name         - [string]              (1x1) the name of the variable
            %variable_idx - [string]              (1x1) the index of the variable
            %
            %OUTPUT
            %
            %none: the index of the requested loading coefficient     
            if not(isempty(obj.X_beta_name))
                index=find(strcmp(obj.X_beta_name{variable_idx},name));
            else
                index=[];
            end
        end   
        
        function index = get_X_z_index(obj,name,variable_idx)
            %DESCRIPTION: returns the index of the variable given its name
            %
            %INPUT
            %obj          - [stem_varset object]  (1x1) the stem_varset object
            %name         - [string]              (1x1) the name of the variable
            %variable_idx - [string]              (1x1) the index of the variable
            %
            %OUTPUT
            %
            %none: the index of the requested loading coefficient    
            if not(isempty(obj.X_z_name))
                index=find(strcmp(obj.X_z_name{variable_idx},name));
            else
                index=[];
            end
        end        
        
        function index = get_X_p_index(obj,name,variable_idx)
            %DESCRIPTION: returns the index of the variable given its name
            %
            %INPUT
            %obj          - [stem_varset object]  (1x1) the stem_varset object
            %name         - [string]              (1x1) the name of the variable
            %variable_idx - [string]              (1x1) the index of the variable
            %
            %OUTPUT
            %
            %none: the index of the requested loading coefficient    
            if not(isempty(obj.X_p_name))
                index=find(strcmp(obj.X_p_name{variable_idx},name));
            else
                index=[];
            end
        end    
        
        function index = get_X_bp_index(obj,name,variable_idx)
            %DESCRIPTION: returns the index of the variable given its name
            %
            %INPUT
            %obj          - [stem_varset object]  (1x1) the stem_varset object
            %name         - [string]              (1x1) the name of the variable
            %variable_idx - [string]              (1x1) the index of the variable
            %
            %OUTPUT
            %
            %none: the index of the requested loading coefficient    
            if not(isempty(obj.X_bp_name))
                index=find(strcmp(obj.X_bp_name{variable_idx},name));
            else
                index=[];
            end
        end         

        %Class set methods
        function nvar = get.nvar(obj)
            nvar=length(obj.Y);
        end
        
        function dim = get.dim(obj)
            dim=zeros(1,length(obj.Y));
            for i=1:length(obj.Y)
                dim(1,i)=size(obj.Y{i},1);
            end
        end
        
        function N = get.N(obj)
            N=0;
            for i=1:length(obj.Y)
                N=N+size(obj.Y{i},1);
            end   
        end
        
        function T = get.T(obj)
            T=size(obj.Y{1},2);
        end
        
        function X_bp_tv = get.X_bp_tv(obj)
            if not(isempty(obj.X_bp))
                if size(obj.X_bp{1},3)==1
                    X_bp_tv=0;
                else
                    X_bp_tv=1;
                end
            else
                X_bp_tv=0;
            end
        end
        
        function X_beta_tv = get.X_beta_tv(obj)
            if not(isempty(obj.X_beta))
                
                if size(obj.X_beta{1},3)==1
                    X_beta_tv=0;
                else
                    X_beta_tv=1;
                end
            else
                X_beta_tv=0;
            end
        end
        
        function X_z_tv = get.X_z_tv(obj) 
            if not(isempty(obj.X_z))
                if size(obj.X_z{1},3)==1
                    X_z_tv=0;
                else
                    X_z_tv=1;
                end
            else
                X_z_tv=0;
            end
        end
        
        function X_p_tv = get.X_p_tv(obj)
            if not(isempty(obj.X_p))
                if size(obj.X_p{1},3)==1
                    X_p_tv=0;
                else
                    X_p_tv=1;
                end
            else
                X_p_tv=0;
            end
        end
        
        function X_h_tv = get.X_h_tv(obj)
            if not(isempty(obj.X_h))
                if size(obj.X_h{1},2)==1
                    X_h_tv=0;
                else
                    X_h_tv=1;
                end
            else
                X_h_tv=0;
            end
        end
        
        function set.Y(obj,Y)
            if not(iscell(Y))
                error('Y must be a cell array');
            end
            obj.Y=Y;
        end

        function set.Y_name(obj,Y_name)
            if not(iscell(Y_name))
                error('Y_name must be a cell array');
            end
            obj.Y_name=Y_name;
        end
        
        function set.X_bp(obj,X_bp)
            if not(iscell(X_bp))
                error('X_bp must be a cell array');
            end
            obj.X_bp=X_bp;
        end  
        
        function set.X_bp_name(obj,X_bp_name)
            if not(iscell(X_bp_name))
                error('X_bp_name must be a cell array');
            end
            obj.X_bp_name=X_bp_name;
        end      
        
        function set.X_beta(obj,X_beta)
            if not(iscell(X_beta))
                error('X_beta must be a cell array');
            end
            obj.X_beta=X_beta;
        end  
        
        function set.X_beta_name(obj,X_beta_name)
            if not(iscell(X_beta_name))
                error('X_beta_name must be a cell array');
            end
            obj.X_beta_name=X_beta_name;
        end      
        
        function set.X_z(obj,X_z)
            if not(iscell(X_z))
                error('X_z must be a cell array');
            end
            obj.X_z=X_z;
        end  
        
        function set.X_z_name(obj,X_z_name)
            if not(iscell(X_z_name))
                error('X_z_name must be a cell array');
            end           
            obj.X_z_name=X_z_name;
        end      
        
        function set.X_p(obj,X_p)
            if not(iscell(X_p))
                error('X_p must be a cell array');
            end
            obj.X_p=X_p;
        end  
        
        function set.X_p_name(obj,X_p_name)
            if not(iscell(X_p_name))
                error('X_p_name must be a cell array');
            end
            obj.X_p_name=X_p_name;
        end          
        
        function set.X_h(obj,X_h)
            if not(iscell(X_h))
                error('X_h must be a cell array');
            end
            obj.X_h=X_h;
        end
        
        function set.X_h_name(obj,X_h_name)
            if not(ischar(X_h_name))
                error('X_h_name must be a string');
            end
            obj.X_h_name=X_h_name;
        end 
    end
    
end