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

classdef stem_misc
    
    methods (Static)
        
        function R = autocorr(y,nlag,flag_plot)
            %DESCRIPTION: evaluate and plot the temporal autocorrelation of the time series y
            %
            %INPUT
            %
            %y                  - [double]      (Nx1)  the time series
            %nlag               - [integer]     (1x1)  number of temporal lags to evaluate
            %flag_plot          - [boolean]     (1x1)  1: plot the result; 0: no plot
            %pixel_side_h=[];   - [double]      (1x1)          the height of the pixel (in the same unit of measure of the unit property)
            %
            %OUTPUT
            %
            %R                  - [double]      (1x1) the autocorrelation evaluation for nlag lags
            if nargin < 2
                nlag = 24;
            end
            if nargin < 3
                flag_plot=1;
            end
            if size(y,2) > 1
                error('Multivariate time series not supported');
            else
                L=isnan(y);
                if not(isempty(L))
                    warning('Missing data removed from the vector');
                    y(L)=[];
                end
                n = length(y);
                mnlag = min(n-1,nlag);
                r = covf( (y - mean(y)), mnlag+1 );
                r = r' / r(1);
                sigma2 = 2 ./ sqrt(n - (0:mnlag));
                sigma3 = 3 ./ sqrt(n - (0:mnlag));
                if flag_plot
                    plot(0:mnlag,r,'-k', 0:mnlag,[sigma2' -sigma2'],'-g', 0:mnlag,[sigma3' -sigma3'],'-r', [0; mnlag], [0; 0]);
                    axis tight
                    title('Autocorrelation (bands at \pm 2 and \pm 3\sigma)' );
                end
            end
            if nargout>0
                R=r;
            end
        end
        
        function res = chol_solve(c,b)
            %DESCRIPTION: linear system solving using Cholesky decomposition - res=a\b
            %
            %INPUT
            %
            %c                  - [double]      (NxN)  the Cholesky decomposition of the s.s.d. matrix a
            %b                  - [double]      (NxM)  the matri b
            %
            %OUTPUT
            %
            %res                - [double]      (NxM)  the result of the linear system a\b
            %solve inv(a)*b where c is chol(a)
            res=c\(c'\b);
        end
        
        function corr = correlation_function(theta,DistMat,correlation_type,idx_r,idx_c)
            %DESCRIPTION: evaluation of the spatial correlation function
            %
            %INPUT
            %
            %theta              - [double]      (1x1)|(2x1) the parameter scalar or vector of the spatial correlation function
            %DistMat            - [double]      (NxN)|{2}(NxN) the distance matrix or matrices
            %correlation_type   - [string]      (1x1) 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2; 'expsphere': anisotropic correlation function on the spehere
            %idx_r              - [integer>0]   (Nx1) the row indexes of the submatrix
            %idx_c              - [integer>0]   (Nx1) the column indexes of the submatrix
            %
            %OUTPUT
            %
            %corr               - [double]      (NxN)   spatial correlation matrix
            if sum(strcmp(correlation_type,{'exponential','matern32','matern52','expsphere'}))==0
                error('The correlation type must be either ''exponential'', ''matern32'' or ''matern52'' or ''expsphere''');
            end
            
            if strcmp(correlation_type,'expsphere')
                if not(iscell(DistMat))
                    error('DistMat must be a 2x1 cell array and each cell a NxN distance matrix');
                else
                    if not(length(DistMat)==2)
                        error('DistMat must be a 2x1 cell array and each cell a NxN distance matrix');
                    end
                end
                if not(size(theta,1)==2)
                    error('theta must be 2x1 or 2xp when type is ''expsphere''');
                end
            else
                if not(length(theta)==1)
                    warning('theta must be a scalar');
                end
            end
            
            if nargin<4
                if not(strcmp(correlation_type,'expsphere'))
                    idx_r=1:size(DistMat,1);
                else
                    idx_r=1:size(DistMat{1},1);
                end
            end
            if nargin<5
                if not(strcmp(correlation_type,'expsphere'))
                    idx_c=1:size(DistMat,2);
                else
                    idx_c=1:size(DistMat{1},2);
                end
            end
            
            if not(strcmp(correlation_type,'expsphere'))
                DistMat=DistMat(idx_r,idx_c);
            else
                DistMat{1}=DistMat{1}(idx_r,idx_c);
                DistMat{2}=DistMat{2}(idx_r,idx_c);
            end
           
            if not(issparse(DistMat))
                if strcmp(correlation_type,'exponential')
                    corr=exp(-DistMat/theta);
                end
                if strcmp(correlation_type,'matern32')
                    corr=(1+sqrt(3)*DistMat/theta).*exp(-sqrt(3)*DistMat/theta);
                end
                if strcmp(correlation_type,'matern52')
                    corr=(1+sqrt(5)*DistMat/theta+5*DistMat.^2/(3*theta^2)).*exp(-sqrt(5)*DistMat/theta);
                end
                if strcmp(correlation_type,'expsphere')
                    corr=exp(-DistMat{1}/theta(1)).*exp(-DistMat{2}/theta(2));
                end
            else
                %note that the diagonal of the variance-covariance matrices and
                %possibly some elements of the cross-covariance matrices are equal
                %to eps instead of zero. However, exp(eps)=1.
                idx=find(DistMat);
                if strcmp(correlation_type,'exponential')
                    correlation=exp(-DistMat(idx)/theta);
                end
                if strcmp(correlation_type,'matern32')
                    correlation=(1+sqrt(3)*DistMat(idx)/theta).*exp(-sqrt(3)*DistMat(idx)/theta);
                end
                if strcmp(correlation_type,'matern52')
                    correlation=(1+sqrt(5)*DistMat(idx)/theta+5*DistMat(idx).^2/(3*theta^2)).*exp(-sqrt(5)*DistMat(idx)/theta);
                end
                if strcmp(correlation_type,'expsphere')
                    correlation=exp(-DistMat{1}(idx)/theta(1)).*exp(-DistMat{2}(idx)/theta(2));
                end
                [I,J]=ind2sub(size(DistMat),idx);
                corr.I=I;
                corr.J=J;
                corr.correlation=correlation;
            end
        end
        
        function res = D_apply(a,d,type)
            %DESCRIPTION: this method is used in EM estimation to avoid full matrix multiplication
            %
            %INPUT
            %
            %a              - [double]      (Nx1|N) vector or matrix
            %d              - [double]      (Mx1)   the vector to pre and/or post multiply. The vector is the diagonal of diagonal matrix
            %type           - [string]      (1x1)   'l': left pre multiplication; 'r': right post multiplication; 'b' pre and post multiplication
            %
            %OUTPUT
            %
            %res            - [double]      (Nx1|N) the result of the multiplication
            
            if nargin<3
                error('All the input arguments must be provided');
            end
            if not(isempty(d))
                if not(sum(abs(d))==0||sum(abs(a(:)))==0)
                    if isvector(a)
                        if strcmp(type,'l')
                            %a is a column vector
                            r=length(d)-length(a);
                            if r==0
                                res=d.*a;
                            else
                                if r>0
                                    if sum(abs(d(end-r+1:end)))==0
                                        res=[d(1:length(a)).*a;zeros(r,1)];
                                    else
                                        error('The elements of d exceeding the dimension of a must be zeros');
                                    end
                                else
                                    res=d.*a;
                                    %error('The vector d cannot be smaller than the vector a');
                                end
                            end
                        end
                        if strcmp(type,'r')
                            %a is a row vector
                            r=length(d)-length(a);
                            if r==0
                                res=a.*d';
                            else
                                if r>0
                                    if sum(abs(d(end-r+1:end)))==0
                                        res=[a.*d(1:length(a))',zeros(1,r)];
                                    else
                                        error('The elements of d exceeding the dimension of a must be zeros');
                                    end
                                else
                                    error('The vector d cannot be smaller than the vector a');
                                end
                            end
                        end
                        if strcmp(type,'b')
                            if not(iscolumn(a))
                                error('type ''b'' D_apply can be used only with ''a'' as a column vector')
                            end
                            r=length(d)-length(a);
                            if r==0
                                res=(d.^2).*a;
                            else
                                if r>0
                                    if sum(abs(d(end-r+1:end)))==0
                                        res=[(d(1:length(a)).^2).*a;zeros(r,1)];
                                    else
                                        error('The elements of d exceeding the dimension of a must be zeros');
                                    end
                                else
                                    error('The vector d cannot be smaller than the vector a');
                                end
                            end
                        end
                    end
                    if min(size(a))>1 %is a matrix
                        if strcmp(type,'l')
                            r=length(d)-size(a,1);
                            if r==0
                                res=sparse(diag(d))*a;
                            else
                                if r>0
                                    if sum(abs(d(end-r+1:end)))==0
                                        res=sparse(diag(d(1:size(a,1))))*a;
                                        res=[res;zeros(r,size(a,2))];
                                    else
                                        res=a.*repmat(d',[size(a,1),1]);
                                    end
                                else
                                    res=a.*repmat(d',[size(a,1),1]);
                                end
                            end
                        else
                            if strcmp(type,'r')
                                r=length(d)-size(a,2);
                                if r==0
                                    res=a*sparse(diag(d));
                                else
                                    if r>0
                                        if sum(abs(d(end-r+1:end)))==0
                                            res=a*sparse(diag(d(1:size(a,2))));
                                            res=[res,zeros(size(a,1),r)];
                                        else
                                            error('The element of d exceeding the dimension of a must be zeros');
                                        end
                                    else
                                        error('The vector d cannot be smaller than the second dimension of the matrix a');
                                    end
                                    
                                end
                            else
                                r=length(d)-size(a,1);
                                if r==0
                                    if not(issparse(a))
                                        res=(d*d').*a;
                                    else
                                        I=1:length(d);
                                        D=sparse(I,I,d);
                                        res=D*a*D;
                                    end
                                else
                                    if r>0
                                        if sum(abs(d(end-r+1:end)))==0
                                            if not(issparse(a))
                                                res=[(d(1:size(a,1))*d(1:size(a,1))').*a, zeros(size(a,1),r); zeros(r,size(a,2)) zeros(r,r)];
                                            else
                                                I=1:length(d(1:size(a,1)));
                                                D=sparse(I,I,d(1:size(a,1)));
                                                res=D*a*D;
                                                L=find(res);
                                                [I,J] = ind2sub(size(res),L);
                                                res=sparse(I,J,full(res(L)),size(res,1)+r,size(res,2)+r);
                                            end
                                        else
                                            error('The element of d exceeding the dimension of a must be zeros');
                                        end
                                    else
                                        error('The vector d cannot be smaller than the dimension of the matrix a');
                                    end
                                end
                            end
                        end
                    end
                else
                    %return a matrix of all zeros
                    if issparse(a)
                        if strcmp(type,'b')
                            res=sparse(length(d),length(d));
                        end
                        if strcmp(type,'l')
                            r=length(d)-size(a,1);
                            if (r==0)
                                res=sparse(length(d),size(a,2));
                            else
                                res=sparse(size(a,1),size(a,2));
                            end
                        end
                        if strcmp(type,'r')
                            res=sparse(size(a,1),length(d));
                        end
                    else
                        if strcmp(type,'b')
                            res=zeros(length(d),length(d));
                        end
                        if strcmp(type,'l')
                            r=length(d)-size(a,1);
                            if r==0
                                res=zeros(length(d),size(a,2));
                            else
                                res=zeros(size(a,1),size(a,2));
                            end
                        end
                        if strcmp(type,'r')
                            res=sparse(size(a,1),length(d));
                        end
                    end
                end
            else
                %nothing to do
                res=a;
            end
        end
        
        function str = decode_time(s)
            %DESCRIPTION: convert an interval of time in seconds into a string with hours, minutes and seconds
            %
            %INPUT
            %
            %s              - [double]      (1x1)   the interval of time in seconds
            %
            %OUTPUT
            %
            %str            - [string]      (1x1)   the output string               
            
            h=floor(s/3600);
            s=s-h*3600;
            m=floor(s/60);
            s=s-m*60;
            if h>0
                str_h=[num2str(h),'h '];
            else
                str_h=[];
            end
            if m>0
                str_m=[num2str(m),'m '];
            else
                str_m=[];
            end
            str_s=[num2str(s,'%5.2f'),'s'];
            str=[str_h,str_m,str_s];
        end
        
        function disp_star(string)
            %DESCRIPTION: embed and display a string into a frame of *
            %
            %INPUT
            %
            %string           - [string]      (1x1)   the string to display
            %
            %OUTPUT
            %
            %none: the string is displayed              
                        
            l=length(string);
            s=[];
            for i=1:l+8
                s=cat(2,s,'*');
            end
            disp(' ');
            disp(s);
            disp(['*   ',string,'   *']);
            disp(s);
        end
        
        function mat = from_vector_to_symmetric_matrix(vec)
            %DESCRIPTION: a correlation matrix is built from the vector vec
            %
            %INPUT
            %
            %vec           - [double]      (dx1)    the vector of the upper extra-diagonal elements
            %
            %OUTPUT
            %
            %mat:          - [double]      (qxq)    the correlation matrix           
            d=(1+sqrt(1+8*length(vec)))/2;
            mat=eye(d);
            counter=1;
            for i=1:d-1
                for j=i+1:d
                    mat(i,j)=vec(counter);
                    mat(j,i)=mat(i,j);
                    counter=counter+1;
                end
            end
            
        end
        
        function vec = from_upper_triangular_to_vector(mat)
            %DESCRIPTION: return the vector of the upper triangular part of a matrix
            %
            %INPUT
            %mat          - [double]      (dxd) the input matrix
            %
            %OUTPUT
            %vec:         - [double]      (d*(d-1)/2x1) the vector of the upper extra-diagonal elements  
            
            d=size(mat,1);
            vec=zeros(d*(d-1)/2,1);
            counter=1;
            for i=1:d-1
                for j=i+1:d
                    vec(counter)=mat(i,j);
                    counter=counter+1;
                end
            end
        end        
        
        function xls_coordinates = rc2xls(row,col)
            if col>26*27
                error('Not supported');
            end
            first_letter=floor((col-1)/26);
            if first_letter>0
                first_letter=char(first_letter+64);
            else
                first_letter='';
            end
            second_letter=char(mod((col-1),26)+1+64);
            xls_coordinates=[first_letter,second_letter,num2str(row)];
        end
        
        function [B,block_i,block_j] = get_block(dim_r,i,dim_c,j,A)
            %DESCRIPTION: returns the block of a block matrix
            %
            %INPUT
            %
            %dim_r           - [integer]      (cx1)   the number of rows in each sub-block
            %i               - [integer]      (1x1)   the row index of the sub-block to extract
            %dim_c           - [integer]      (dx1)   the number of columns in each sub-block
            %j               - [integer]      (1x1)   the column index of the sub-block to extract
            %A               - [double]       (NxN)   the block matrix
            %
            %OUTPUT
            %
            %B               - [double]       (GxH)   the matrix sub-block
            %block_i         - [integer]      (Gx1)   the indices of the extracted rows
            %block_j         - [integer]      (Hx1)   the indices of the extracted columns
            
            rr=cumsum(dim_r);
            block_i = rr(i)-dim_r(i)+1:rr(i);
            if nargin<=2
                % rende solo block_i
                B=block_i;
                block_i=[];
                block_j=[];
                return
            end
            
            cc=cumsum(dim_c);
            block_j = cc(j)-dim_c(j)+1:cc(j);
            if nargin<=4
                B=block_i;
                block_i=block_j;
                block_j=[];
                return
            end
            B=A(block_i,block_j);
        end
        
        function color = get_rainbow_color(value,min_value,max_value)
            %DESCRIPTION: given a value within a range, the respective color in the rainbow colormap is returned
            %
            %INPUT
            %
            %value          - [double]       (1x1)   the value within the range
            %min_value      - [double]       (1x1)   the lower limit of the range
            %max_value      - [double]       (1x1)   the upper limit of the range
            %
            %OUTPUT
            %
            %color          - [double >=0 and <=1] (3x1)    the RGB color 
            
            value=(value-min_value)/(max_value-min_value);
            if value<0.25
                color(1)=0;
                color(2)=value*4;
                color(3)=1;
            else
                if value<0.5
                    color(1)=0;
                    color(2)=1;
                    color(3)=1-(value-0.25)*4;
                else
                    if value<0.75
                        color(1)=(value-0.5)*4;
                        color(2)=1;
                        color(3)=0;
                    else
                        color(1)=1;
                        color(2)=1-(value-0.75)*4;
                        color(3)=0;
                    end
                end
            end
        end
            
        function result = isdiagonal(a)
            %DESCRIPTION: test if a matrix is diagonal
            %
            %INPUT
            %
            %a          - [double]       (NxN)   the input matrix
            %
            %OUTPUT
            %
            %result     - [boolean]      (1x1)   1: the matrix is diagonal; 0: otherwise
            [i,j] = find(a);
            if ~isempty(i)
                result = all(i == j);
            else
                result = true;
            end
        end
        
        function res = M_apply(a,M,type)
            %DESCRIPTION: this method is used in EM estimation to avoid full matrix multiplication
            %
            %INPUT
            %
            %a              - [double]      (Nx1|N) vector or matrix
            %M              - [integer >0]  (N_gx1) vector of indices
            %type           - [string]      (1x1)   'l': left pre multiplication; 'r': right post multiplication; 'b' pre and post multiplication
            %
            %OUTPUT
            %
            %res            - [double]      (Nx1|N) the result of the multiplication            

            if nargin<3
                error('All the input arguments must be provided');
            end
            if isvector(a)
                if strcmp(type,'l')
                    res=[a(M);a];
                end
                if strcmp(type,'r')
                    res=[a(M) a];
                end
            end
            if ismatrix(a)
                if strcmp(type,'l')
                    if not(issparse(a))
                        res=[a(M,:);a];
                    else
                        res2=a(M,:);
                        L2=find(res2);
                        [I2,J2] = ind2sub(size(res2),L2);
                        res2=full(res2(L2));
                        
                        res3=a;
                        L3=find(res3);
                        [I3,J3] = ind2sub(size(res3),L3);
                        I3=I3+length(M);
                        res3=full(res3(L3));
                        
                        res=[res2;res3];
                        I=[I2;I3];
                        J=[J2;J3];
                        res=sparse(I,J,res);
                    end
                else
                    if strcmp(type,'r')
                        if not(issparse(a))
                            res=[a(:,M) a];
                        else
                            res2=a(:,M);
                            L2=find(res2);
                            [I2,J2] = ind2sub(size(res2),L2);
                            res2=full(res2(L2));
                            
                            res3=a;
                            L3=find(res3);
                            [I3,J3] = ind2sub(size(res3),L3);
                            J3=J3+length(M);
                            res3=full(res3(L3));
                            
                            res=[res2;res3];
                            I=[I2;I3];
                            J=[J2;J3];
                            res=sparse(I,J,res);
                        end
                    else
                        if not(issparse(a))
                            res=[a(M,M) a(M,:); a(:,M) a];
                        else
                            res1=a(M,M);
                            L1=find(res1);
                            [I1,J1] = ind2sub(size(res1),L1);
                            res1=full(res1(L1));
                            
                            res2=a(M,:);
                            L2=find(res2);
                            [I2,J2] = ind2sub(size(res2),L2);
                            J2=J2+length(M);
                            res2=full(res2(L2));
                            
                            res3=a;
                            L3=find(res3);
                            [I3,J3] = ind2sub(size(res3),L3);
                            I3=I3+length(M);
                            J3=J3+length(M);
                            res3=full(res3(L3));
                            
                            res=[res1;res2;res2;res3];
                            I=[I1;I2;J2;I3];
                            J=[J1;J2;I2;J3];
                            res=sparse(I,J,res);
                        end
                    end
                end
            end
        end
        
        function h_fig = plot_map(lat,lon,data,shape,fig_title,fig_xlabel,fig_ylabel)
            %DESCRIPTION: plot geolocated data over a map (if given)
            %
            %INPUT
            %
            %lat            - [double]      (Nx1) vector of latitude
            %lon            - [double]      (Nx1) vector of longitude
            %data           - [double]      (Nx1) the data to plot
            %<shape>        - [shape file]  (1x1) the shape file of the map
            %<fig_title>    - [string]      (1x1) (default: []) the title of the figure
            %<fig_xlabel>   - [string]      (1x1) (default: []) the xlabel of the figure
            %<fig_ylabel>   - [string]      (1x1) (default: []) the ylabel of the figure
            %
            %OUTPUT
            %
            %none: the data are plotted
            
            if nargin<3
                error('lat, lon and data must be provided');
            end
            if nargin<4
                shape=[];
            end
            if nargin<5
                fig_title=[];
            end
            if nargin<6
                fig_xlabel=[];
            end
            if nargin<7
                fig_ylabel=[];
            end
            if not(isvector(lat))&&not(isvector(lon))&&not(isvector(data))
                h=figure;
                hold on
                if not(isempty(shape))
                    mapshow(shape);
                end
                h2 = mapshow(lon,lat,data,'DisplayType','texture');
                set(h2,'FaceColor','flat');
                axis equal
                xlim([min(lon(:)),max(lon(:))]);
                ylim([min(lat(:)),max(lat(:))]);
                colormap summer;
                colorbar;                
                title(fig_title,'FontSize',16);
                xlabel(fig_xlabel,'FontSize',16);
                ylabel(fig_ylabel,'FontSize',16);
                grid on;
                box on;
                set(gca,'FontSize',16);
                set(gcf, 'renderer', 'zbuffer');
            else
                if isvector(lat)&&isvector(lon)&&isvector(data)
                    h=figure;
                    if not(isempty(shape))
                        latmin=min(lat);
                        latmax=max(lat);
                        lonmin=min(lon);
                        lonmax=max(lon);
                        rlat=latmax-latmin;
                        rlon=lonmax-lonmin;
                        mapshow(shape);
                        axis equal
                        xlim([lonmin-rlon*0.1,lonmax+rlon*0.1]);
                        ylim([latmin-rlat*0.1,latmax+rlat*0.1]);
                        hold on
                    end
                    minval=nanmin(data);
                    maxval=nanmax(data);
                    for i=1:length(lat)
                        if not(isnan(data(i)))
                            color = stem_misc.get_rainbow_color(data(i),minval,maxval);
                            mapshow(lon(i),lat(i),'DisplayType','point','MarkerFaceColor',color, 'MarkerEdgeColor','k','Marker','o','MarkerSize',8);
                            hold on
                        end
                    end
                    title(fig_title,'FontSize',16);
                    xlabel(fig_xlabel,'FontSize',16);
                    ylabel(fig_ylabel,'FontSize',16);
                    grid on;
                    box on;
                    set(gca,'FontSize',16);
                    set(gcf, 'renderer', 'zbuffer');
                else
                    error('lat, lon and data must be either all matrices or all vectors');
                end
            end
            if nargout>0
                h_fig=h;
            end
        end
        
        function m = sparseif(m,density)
            %DESCRIPTION: return the sparse matrix of the matrix m if the zero density is higher than density otherwise the full matrix
            %
            %INPUT
            %
            %m                 - [double]      (dxd) the input matrix
            %density           - [double]      (1x1) the zero density of matrix m expressed as percentage
            %
            %OUTPUT
            %
            %m                 - [double]      (dxd) the output matrix      
            if stem_misc.zero_density(m)>density
                m=sparse(m);
            else
                m=full(m);
            end
        end
        
        function v = triuv(mat)
            %DESCRIPTION: extract the diagonal and the upper triangular part of a matrix into a vector
            %
            %INPUT
            %
            %mat            - [double]      (dxd) the input matrix
            %
            %OUTPUT
            %
            %v              - [double]      (d*(d+1)/2x1) the output vector
            
            d=size(mat,1);
            v=zeros(d*(d+1)/2,1);
            l=1;
            for j=1:d
                for i=j:d
                    v(l)=mat(j,i);
                    l=l+1;
                end
            end
        end
        
        function weights = wendland(DistMat,gamma,correlation_type,flag_reshape,idx_r,idx_c)
            %DESCRIPTION: compute the tapering weights using the Wendland function
            %
            %INPUT
            %
            %DistMat           - [double]      (NxN)    the distance matrix
            %gamma             - [double >0]   (1x1)    the gamma parameter of the wendland function
            %correlation_type  - [string]      (1x1)    'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2; 'expsphere': anisotropic correlation function on the spehere
            %flag_reshape      - [boolean]     (1x1)    1:the weights are reshaped to have the same dimension of the distance matrix; 0: the weights are given only for the non-zero elements
            %idx_r             - [integer>0]   (Nx1)    the row indexes of the submatrix
            %idx_c             - [integer>0]   (Nx1)    the column indexes of the submatrix
            %
            %OUTPUT
            %
            %weights:          - [double]      (NxN | dNx1) the tapering weights. The dimension of the output depends on the flag_reshape value
            
            
            if nargin<3
                error('The first 3 inputs arguments must be provided');
            end
            if nargin<4
                flag_reshape=0;
            end
            if gamma<=0
                error('The gamma tapering parameter must be > 0');
            end
            
            if nargin<5
                if not(strcmp(correlation_type,'expsphere'))
                    idx_r=1:size(DistMat,1);
                else
                    idx_r=1:size(DistMat{1},1);
                end
            end
            if nargin<6
                if not(strcmp(correlation_type,'expsphere'))
                    idx_c=1:size(DistMat,2);
                else
                    idx_c=1:size(DistMat{1},2);
                end
            end
            
            if not(strcmp(correlation_type,'expsphere'))
                DistMat=DistMat(idx_r,idx_c);
            else
                DistMat=DistMat{1}(idx_r,idx_c); %only the geodetic distance is retained
            end
            
            %find the indices of the non-zero elements of DistMat;
            idx=find(DistMat);
            %evaluate the wendland function only on the non-zero elements
            c=DistMat(idx)./gamma;
            weights=((1-c).^4).*(4*c+1); %or %weights=((1-c).^6).*(35*c.^2+18*c+3);
            if flag_reshape
                [I,J]=ind2sub(size(DistMat),idx);
                weights=sparse(I,J,weights,size(DistMat,1),size(DistMat,2));
            end
        end
        
        function compare(ref_par,par,name)
            if not(isempty(ref_par))
                size_ref_par=size(ref_par);
                size_par=size(par);
                
                if not(length(size_ref_par)==length(size_par))
                    error([name,' must be ',num2str(size_ref_par(1)),'x',num2str(size_ref_par(2))]);
                else
                    if sum(abs(size_ref_par-size_par))>0
                        if length(size_ref_par)==2
                            warning([name,' must be ',num2str(size_ref_par(1)),'x',num2str(size_ref_par(2))]);
                        else
                            error([name,' must be ',num2str(size_ref_par(1)),'x',num2str(size_ref_par(2)),'x',num2str(size_ref_par(3))]);
                        end
                    end
                end
                
            end
        end
        
        function name = varname(var)
            name = inputname(1);
        end
        
        function v = vec(mat)
            %DESCRIPTION: vectorize a matrix
            %
            %INPUT
            %
            %mat             - [double]       (Nx...xM)   the input matrix
            %
            %OUTPUT
            %
            %v              - [double]        (N*...*Mx1) the vectorize matrix
            v=mat(:);
        end
        
        function density = zero_density(matrix)
            %DESCRIPTION: returns the density of zero elements within the matrix
            %
            %INPUT
            %
            %matrix          - [double]       (NxN)   the input matrix
            %
            %OUTPUT
            %
            %density         - [double]       (1x1)   the percentage of zero elements            
            density=100-nnz(matrix)/(size(matrix,1)*size(matrix,2))*100;
        end
        
        %Yaqiong
        function f = update_coe_sigma_eps(x,i,others,Omega,Basis,flag_logsigma)
            %DESCRIPTION: update sigma_eps when modeltype is f-HDGM and we
            %set splines to estimate sigma_eps, the function will be called
            %at stem_EM.M_setp when updating basis coefficients of sigma_eps
            %
            %INPUT
            %
            %Omega          - [double]       (N*Tx1)   the input matrix
            %Basis          - [double]       (N*Txk)   the input matrix
            %coe_sigma_eps  - [double]        (kx1) the vector
            %
            %OUTPUT
            %
            %n = size(Omega,1)/T;
            others(i) = x;
            coe_spline_sigma_eps = others;
            f = 0;
            for t=1:size(Omega,1)
                if flag_logsigma==1
                    sigma_t = exp(Basis(t,:)*coe_spline_sigma_eps);
                else
                    sigma_t = (Basis(t,:)*coe_spline_sigma_eps).^2;
                end
                f = f + log(sigma_t)+(sigma_t)^(-1)*Omega(t);
            end  
        end
        
        function [obj_stem_varset,obj_stem_gridlist_p]=table_to_varset(DataTable)

            warning('the function is only avaliable ');
            LatLon=unique([DataTable.Lat,DataTable.Lon],'rows');
            N=size(LatLon,1);
            T=max(DataTable.Time_step);
            k=strfind(DataTable.Properties.VariableNames,'X_f');
            for varidx=1:size(DataTable,2)
                if k{varidx}
                    Max_q=length(DataTable{1,varidx}{1});
                    for i=2:size(DataTable,1)
                        Max_q=max(Max_q,length(DataTable{i,varidx}{1}));
                    end
                end
            end
            

            Y=cell(Max_q);
            Y_name=cell(Max_q);
            X_beta=cell(Max_q);
            X_beta_name=cell(Max_q);
            X_f=cell(Max_q);
            X_f_name=cell(Max_q);
            for i=1:Max_q
                %response variable Y
                k=strfind(DataTable.Properties.VariableNames,'Y_');
                for varidx=1:size(DataTable,2)
                    if k{varidx}
                        for s=1:N
                            for t=1:T
                                tmp=table2array(DataTable(DataTable.Time_step==t&DataTable.Lat==LatLon(s,1)&DataTable.Lon==LatLon(s,2),varidx));
                                Y{i}(s,t)=tmp{1}(i);
                            end
                        end 
                        tmp=strsplit(DataTable.Properties.VariableNames{varidx},'_');
                        Y_name{i}=tmp(end);
                    end
                end

                %X_f
                k=strfind(DataTable.Properties.VariableNames,'X_f');
                for varidx=1:size(DataTable,2)
                    if k{varidx}
                        for s=1:N
                            for t=1:T
                                tmp=table2array(DataTable(DataTable.Time_step==t&DataTable.Lat==LatLon(s,1)&DataTable.Lon==LatLon(s,2),varidx));
                                X_f{i}(s,t)=tmp{1}(i);
                            end
                        end 
                        tmp=strsplit(DataTable.Properties.VariableNames{varidx},'_');
                        X_f_name{i}=tmp(end);
                    end
                end

                %X_beta
                k=strfind(DataTable.Properties.VariableNames,'X_beta');
                counter=1;
                X_beta_name_tmp=[];
                for varidx=1:size(DataTable,2)
                    if k{varidx}
                        for s=1:N
                            for t=1:T
                                tmp=table2array(DataTable(DataTable.Time_step==t&DataTable.Lat==LatLon(s,1)&DataTable.Lon==LatLon(s,2),varidx));
                                X_beta{i}(s,counter,T)=tmp{1}(i);
                            end
                        end 
                        counter=counter+1;
                        tmp=strsplit(DataTable.Properties.VariableNames{varidx},'_');
                        X_beta_name_tmp=cat(2,X_beta_name_tmp,tmp(end));
                    end
                end
                X_beta_name{i}=X_beta_name_tmp;
                clear X_beta_name_tmp
            end
            X_bp=[];
            X_bp_name=[];
            X_z=[];
            X_z_name=[];
            X_p=[];
            X_p_name=[];
            obj_stem_varset=stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name,X_f,X_f_name);

            obj_stem_gridlist_p=stem_gridlist();
            obj_stem_grid=stem_grid(LatLon,'deg','sparse','point');
            obj_stem_grid=obj_stem_grid.sorted_by_lat;
            for i=1:q  
                obj_stem_gridlist_p.add(obj_stem_grid);
            end
        end
      

    end
end