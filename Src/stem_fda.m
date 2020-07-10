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
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_fda < handle
    
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
    %s - dimension of spline knots
   
    properties (SetAccess=private)
        spline_type=[]              %[string]       (1x1) the type of basis, either 'Bspline' or 'Fourier'
        spline_range=[];            %[double]       (2x1) the range of the basis function
        
        spline_order_beta=[];       %[integer >0]   (1x1) the order of the spline for beta 
        spline_knots_beta=[];       %[double]       (sx1) the spline knots for beta 
        spline_basis_beta=[];       %[basis obj]    (1x1) the spline basis for beta 
        spline_nbasis_beta =[]      %[integer >0]   (1x1) the number of basis for beta 
 
        spline_order_z=[];          %[integer >0]   (1x1) the order of the spline 
        spline_knots_z=[];          %[double]       (sx1) the spline knots 
        spline_basis_z=[];          %[basis obj]    (1x1) the spline basis 
        spline_nbasis_z =[]         %[integer >0]   (1x1) the number of basis 
        
        spline_order_sigma=[];      %[integer >0]   (1x1) the order of the spline for sigma_eps2 
        spline_knots_sigma=[];      %[double]       (sx1) the spline knots for sigma_eps2 
        spline_basis_sigma=[];      %[basis obj]    (1x1) the spline basis for sigma_eps2 
        spline_nbasis_sigma =[]     %[integer >0]   (1x1) the number of basis for sigma_eps2 
        
        flag_beta_spline = 0;       %[interger]     (1x1) the indicator of estimating beta(h)
        flag_sigma_eps_spline = 0;  %[interger]     (1x1) the indicator of estimating sigma_eps2 (h)
    end
    
    properties (Hidden=true)
        flag_logsigma = 1;         %[interger] (1x1) the indicator of is on logsigma or not, with 1 indicating using the log and can not be set by users.
        flag_sqrsigma = 0;         %[interger] (1x1) the indicator of is on sqrtsigma or not, with 1 indicating using the sqrt and can not be set by users.
    end
    
    methods
        
        function obj = stem_fda(struct_fda)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %struct_fda = []           - [struct]            (1x1) structure containing spline information
            %
            %The struct including the following
            %spline_type=[]            - [string]            (1x1) the type of basis, either 'Bspline' or 'Fourier'
            %spline_range =[]          - [double]            (2x1) the range of the basis function
            %
            %when type is Fourier
            %spline_nbasis_z =[]       - [integer >0]        (1x1) the number of basis 
            %spline_nbasis_beta =[]    - [integer >0]        (1x1) the number of basis for beta
            %spline_nbasis_sigma =[]   - [integer >0]        (1x1) the number of basis for sigma_eps
            %
            %when type is Bspline
            %spline_order_z=[]         - [integer >0]        (1x1) the order of the spline
            %spline_knots_z=[]         - [double]            (sx1) the spline knots
            %spline_order_beta=[]      - [integer >0]        (1x1) the order of the spline
            %spline_knots_beta=[]      - [double]            (sx1) the spline knots
            %spline_order_sigma=[]     - [integer >0]        (1x1) the order of the spline
            %spline_knots_sigma=[]     - [double]            (sx1) the spline knots
            %
            %OUTPUT
            %obj                       - [stem_fda object]   (1x1)
            
            path = mfilename('fullpath');
            filepath = fileparts(path);
            path = [filepath,'/fda'];
            fda_file = [filepath,'/fda/create_bspline_basis.m'];
            
            if ~exist(path, 'dir') || ~exist(fda_file, 'file')
                disp('fda toolbox not found. Downloading from FDA webpage by Jim Ramsay...');
                mkdir(path);
                options = weboptions('Timeout',Inf);
                websave([filepath,'/fda/fdaM.zip'],'http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/fdaM.zip',options);
                disp('Unzipping fda toolbox...');
                unzip([filepath,'/fda/fdaM.zip'],[filepath,'/fda']);
                addpath(genpath(filepath));
                disp('fda installation completed.');
            end
            
            if length(fieldnames(struct_fda))<3
                error('All the input arguments must be provided');
            end
            if sum(strcmp(struct_fda.spline_type,{'Fourier','Bspline'}))==0
                error('The spline type must be Fourier or Bspline');
            end
            
            if sum(strcmp('spline_nbasis',fieldnames(struct_fda)))
                struct_fda.spline_nbasis_z=struct_fda.spline_nbasis;
                struct_fda.spline_nbasis_beta=struct_fda.spline_nbasis;
                struct_fda.spline_nbasis_sigma=struct_fda.spline_nbasis;
            end
            if sum(strcmp('spline_order',fieldnames(struct_fda)))
                struct_fda.spline_order_z=struct_fda.spline_order;
                struct_fda.spline_order_beta=struct_fda.spline_order;
                struct_fda.spline_order_sigma=struct_fda.spline_order;
            end
            if sum(strcmp('spline_knots',fieldnames(struct_fda)))
                struct_fda.spline_knots_z=struct_fda.spline_knots;
                struct_fda.spline_knots_beta=struct_fda.spline_knots;
                struct_fda.spline_knots_sigma=struct_fda.spline_knots;
            end
            
            obj.spline_type = struct_fda.spline_type;
            if strcmp(struct_fda.spline_type,'Fourier')
                
                if isempty(struct_fda.spline_nbasis_z)||isempty(struct_fda.spline_range)
                    error('All the input arguments must be provided');
                end
                obj.spline_range=struct_fda.spline_range;
                
                if mod(struct_fda.spline_nbasis_z,2)==0
                    warning('The Fourier basis number must be odd, the number is increased by 1')
                end
                obj.spline_basis_z=create_fourier_basis(obj.spline_range, struct_fda.spline_nbasis_z);
                obj.spline_nbasis_z=getnbasis(obj.spline_basis_z);
                
                if not(isempty(struct_fda.spline_nbasis_beta))
                    obj.flag_beta_spline = 1;
                    
                    if mod(struct_fda.spline_nbasis_beta,2)==0
                        warning('The Fourier basis number must be odd, the number is increased by 1.')
                    end
                    obj.spline_basis_beta=create_fourier_basis(obj.spline_range, struct_fda.spline_nbasis_beta);
                    obj.spline_nbasis_beta=getnbasis(obj.spline_basis_beta);
                end
                
                if not(isempty(struct_fda.spline_nbasis_sigma))
                    obj.flag_sigma_eps_spline = 1;
                    
                    if mod(struct_fda.spline_nbasis_sigma,2)==0
                        warning('The Fourier basis number must be odd, the number is increased by 1.')
                    end
                    obj.spline_basis_sigma=create_fourier_basis(obj.spline_range, struct_fda.spline_nbasis_sigma);
                    obj.spline_nbasis_sigma=getnbasis(obj.spline_basis_sigma);
                end
                
            else
                if isempty(struct_fda.spline_knots_z)||isempty(struct_fda.spline_range)
                    error('Spline knots for z must be provided');
                end
                if max(struct_fda.spline_knots_z)>struct_fda.spline_range(2)
                    error('The highest element of spline_knots cannot be higher than spline_range(2)');
                end
                if min(struct_fda.spline_knots_z)<struct_fda.spline_range(1)
                    error('The lowest element of spline_knots cannot be lower than spline_range(1)');
                end
                
                obj.spline_range = struct_fda.spline_range;
                obj.spline_order_z=struct_fda.spline_order_z;
                obj.spline_knots_z=struct_fda.spline_knots_z;

                norder=obj.spline_order_z+1;
                nbasis=length(obj.spline_knots_z)+norder-2;
                obj.spline_basis_z=create_bspline_basis(obj.spline_range, nbasis, norder, obj.spline_knots_z);
                obj.spline_nbasis_z=getnbasis(obj.spline_basis_z);
                
                if not(isempty(struct_fda.spline_order_beta))&&isempty(struct_fda.spline_knots_beta)
                    error('Spline knots for beta must be provided');
                end
                if (isempty(struct_fda.spline_order_beta))&&not(isempty(struct_fda.spline_knots_beta))
                    error('Spline order for beta must be provided');
                end
                if not(isempty(struct_fda.spline_order_sigma))&&isempty(struct_fda.spline_knots_sigma)
                    error('Spline knots for sigma_eps must be provided');
                end
                if (isempty(struct_fda.spline_order_sigma))&&not(isempty(struct_fda.spline_knots_sigma))
                    error('Spline order for sigma_eps must be provided');
                end

                if not(isempty(struct_fda.spline_order_beta))
                    
                    obj.flag_beta_spline = 1;
                    if max(struct_fda.spline_knots_beta)>struct_fda.spline_range(2)
                        error('The highest element of spline_knots cannot be higher than spline_range(2)');
                    end
                    if min(struct_fda.spline_knots_beta)<struct_fda.spline_range(1)
                        error('The lowest element of spline_knots cannot be lower than spline_range(1)');
                    end

                    obj.spline_order_beta = struct_fda.spline_order_beta;
                    obj.spline_knots_beta = struct_fda.spline_knots_beta;
                    
                    norder=obj.spline_order_beta+1;
                    nbasis=length(obj.spline_knots_beta)+norder-2;
                    obj.spline_basis_beta=create_bspline_basis(obj.spline_range, nbasis, norder, obj.spline_knots_beta);
                    obj.spline_nbasis_beta=getnbasis(obj.spline_basis_beta);
                end
                if not(isempty(struct_fda.spline_order_sigma))
                    
                    obj.flag_sigma_eps_spline = 1;
                    if max(struct_fda.spline_knots_sigma)>struct_fda.spline_range(2)
                        error('The highest element of spline_knots cannot be higher than spline_range(2)');
                    end
                    if min(struct_fda.spline_knots_sigma)<struct_fda.spline_range(1)
                        error('The lowest element of spline_knots cannot be lower than spline_range(1)');
                    end
                    
                    obj.spline_order_sigma = struct_fda.spline_order_sigma;
                    obj.spline_knots_sigma = struct_fda.spline_knots_sigma;

                    norder=obj.spline_order_sigma+1;
                    nbasis=length(obj.spline_knots_sigma)+norder-2;
                    obj.spline_basis_sigma=create_bspline_basis(obj.spline_range, nbasis, norder, obj.spline_knots_sigma);
                    obj.spline_nbasis_sigma=getnbasis(obj.spline_basis_sigma);
                end
            end    
        end
             
        function n = get_basis_number(obj)
            n.beta=[];
            n.z=[];
            n.sigma=[];
            
            if not(isempty(obj.spline_basis_beta))
                n.beta=getnbasis(obj.spline_basis_beta);
            end
            if not(isempty(obj.spline_basis_z))
                n.z=getnbasis(obj.spline_basis_z);
            end
            if not(isempty(obj.spline_basis_sigma))
                n.sigma=getnbasis(obj.spline_basis_sigma);
            end
        end
        
        %Class set methods    
        function set.spline_order_z(obj,spline_order)
            if spline_order<1
                error('spline_order must be > 0');
            end
            obj.spline_order_z=spline_order;
        end
        
        function set.spline_range(obj,spline_range)
            if not(length(spline_range)==2)
                error('spline_range must be a 2x1 vector');
            end
            if spline_range(2)<=spline_range(1)
                error('The upper bound of spline_range must be higher than the lower bound');
            end
            obj.spline_range=spline_range;
        end
        
        function set.spline_knots_z(obj,spline_knots)
            if length(spline_knots)<2
                error('spline_knots must be at least 2x1');
            end
            if any(diff(spline_knots)<=0)
                error('spline_knots must be a vector of sorted and non equal elements');
            end
            obj.spline_knots_z=spline_knots;
        end
        
        function set.spline_order_beta(obj,spline_order_beta)
            if spline_order_beta<1
                error('spline_order for beta must be > 0');
            end
            obj.spline_order_beta=spline_order_beta;
        end
        
        function set.spline_order_sigma(obj,spline_order_sigma)
            if spline_order_sigma<1
                error('spline_order for sigma must be > 0');
            end
            obj.spline_order_sigma=spline_order_sigma;
        end
        
        function set.spline_knots_beta(obj,spline_knots_beta)
            if length(spline_knots_beta)<2
                error('spline_knots for beta must be at least 2x1');
            end
            if any(diff(spline_knots_beta)<=0)
                error('spline_knots for beta must be a vector of sorted and non equal elements');
            end
            obj.spline_knots_beta=spline_knots_beta;
        end
        
        function set.spline_knots_sigma(obj,spline_knots_sigma)
            if length(spline_knots_sigma)<2
                error('spline_knots for sigma must be at least 2x1');
            end
            if any(diff(spline_knots_sigma)<=0)
                error('spline_knots for sigma must be a vector of sorted and non equal elements');
            end
            obj.spline_knots_sigma=spline_knots_sigma;
        end
        
    end
end
