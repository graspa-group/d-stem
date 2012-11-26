%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_misc
    
    methods (Static)
        
        function R = autocorr(y,nlag,fl_plot)
            if nargin < 2
                nlag = 24;
            end
            if nargin < 3
                fl_plot=1;
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
                if fl_plot
                    plot(0:mnlag,r,'-k', 0:mnlag,[sigma2' -sigma2'],'-g', 0:mnlag,[sigma3' -sigma3'],'-r', [0; mnlag], [0; 0]);
                    axis tight
                    title('Autocorrelation (bands at \pm 2 and \pm 3\sigma)' );
                end
            end
            if nargout>0
                R=r;
            end
        end
        
        function res = chol_solve(c,b,trim)
            if nargin<3
                trim=0;
            end
            %solve inv(a)*b where c is chol(a)
            if not(issparse(c))||(trim==0)
                res=c\(c'\b);
                %is 19% faster than
                %y=(c'\b)
                %res=c\y
            else
                y=(c'\b);
                [I,J]=find(abs(y)>1e-3);
                L=sub2ind(size(y),I,J);
                elements=y(L);
                y=sparse(I,J,elements,size(y,1),size(y,2));
                res=c\y;
                [I,J]=find(abs(res)>1e-3);
                L=sub2ind(size(res),I,J);
                elements=res(L);
                res=sparse(I,J,elements,size(res,1),size(res,2));
            end
        end
        
        function corr = correlation_function(theta,DistMat,type)
            if strcmp(type,'exponential')
                %theta(1)=theta
                if not(isscalar(theta))
                    error('Theta must be a scalar');
                end
                if not(issparse(DistMat))
                    corr=exp(-DistMat/theta);
                else
                    %note that the diagonal of the variance-covariance matrices and
                    %possibly some elements of the cross-covariance matrices are equal
                    %to eps instead of zero. However, exp(eps)=1.
                    idx=find(DistMat);
                    correlation=exp(-DistMat(idx)/theta);
                    [I,J]=ind2sub(size(DistMat),idx);
                    corr.I=I;
                    corr.J=J;
                    corr.correlation=correlation;
                end
            else
                %thata(1)=alpha
                %theta(2)=nu
                if issparse(DistMat)
                    error('Sparse matrices not yet supported with Matern correlation function');
                end
                if not(size(theta,2)==2)
                    error('Theta must be a 2x1 vector');
                end
                corr=1/(2^(theta(2)-1)*gamma(theta(2)))*((DistMat/theta(1)).^theta(2)).*besselk(theta(2),DistMat/theta(1));
            end
        end
        
        function res = D_apply(a,d,type)
            %res=d*a
            %a must be either a vector or a matrix!
            %d must be a column vector representing the diagonal of a diagonal matrix
            %type must be 'l' for left, 'r' for right or 'b' for both
            if nargin<3
                error('All the input arguments must be provided');
            end
            if not(sum(d-ones(size(d)))==0)
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
                    if issparse(a)
                        nonzero=nnz(a);
                    end
                    if strcmp(type,'l')
                        r=length(d)-size(a,1);
                        if r==0
                            if not(issparse(a))
                                res=a;
                                for i=1:size(a,1)
                                    res(i,:)=res(i,:)*d(i);
                                end
                            else
                                I=1:length(d);
                                D=sparse(I,I,d);
                                res=D*a;
                            end
                        else
                            if r>0
                                if sum(abs(d(end-r+1:end)))==0
                                    if not(issparse(a))
                                        res=a;
                                        for i=1:size(a,1)
                                            res(i,:)=res(i,:)*d(i);
                                        end
                                        res=[res;zeros(r,size(a,2))];
                                    else
                                        I=1:length(d(1:size(a,1)));
                                        D=sparse(I,I,d(1:size(a,1)));
                                        res=D*a;
                                        L=find(res);
                                        [I,J] = ind2sub(size(res),L);
                                        res=sparse(I,J,full(res(L)),size(res,1)+r,size(res,2));
                                    end
                                else
                                    error('The element of d exceeding the dimension of a must be zeros');
                                end
                            else
                                error('The vector d cannot be smaller than the first dimension of the matrix a');
                            end
                        end
                    else
                        if strcmp(type,'r')
                            r=length(d)-size(a,2);
                            if r==0
                                if not(issparse(a))
                                    res=a;
                                    for i=1:size(a,2)
                                        res(:,i)=res(:,i)*d(i);
                                    end
                                else
                                    I=1:length(d);
                                    D=sparse(I,I,d);
                                    res=a*D;
                                end
                            else
                                if r>0
                                    if sum(abs(d(end-r+1:end)))==0
                                        if not(issparse(a))
                                            res=a;
                                            for i=1:size(a,2)
                                                res(:,i)=res(:,i)*d(i);
                                            end
                                            res=[res,zeros(size(a,1),r)];
                                        else
                                            I=1:length(d(1:size(a,1)));
                                            D=sparse(I,I,d(1:size(a,1)));
                                            res=a*D;
                                            L=find(res);
                                            [I,J] = ind2sub(size(res),L);
                                            res=sparse(I,J,full(res(L)),size(res,1),size(res,2)+r);
                                        end
                                        
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
                if not(length(d)==max(size(a)))
                    error('The vector d of all ones must have the same dimension of a');
                end
                res=a;
            end
            
        end
        
        function str = decode_time(s)
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
            l=length(string);
            s=[];
            for i=1:l+8
                s=[s,'*'];
            end
            disp(' ');
            disp(s);
            disp(['*   ',string,'   *']);
            disp(s);
        end
        
        function mat = from_vector_to_symmetric_matrix(vec)
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
        
        function [B,block_i,block_j] = get_block(dim_r,i,dim_c,j,A)
            % [B,block_i,block,j] = get_block(dim_r,i,dim_c,j,A)
            %
            % rende B=A(block_i,block,j)
            % dim_r = struttura a blocchi delle righe
            % dim_c = struttura a blocchi delle colonne
            
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
            [i,j] = find(a);
            if ~isempty(i)
                result = all(i == j);
            else
                result = true;
            end
        end
        
        function res = M_apply(a,M,type)
            %a must be either a vector or a simmetric matrix!
            %M must be a vector of indices
            %type must be 'l' for left, 'r' for right or 'b' for both
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
        
        function plot_map(lat,lon,data,shape)
            if nargin<3
                error('lat, lon and data must be provided');
            end
            if nargin<4
                shape=[];
            end
            if not(isvector(lat))&&not(isvector(lon))&&not(isvector(data))
                figure
            else
                if isvector(lat)&&isvector(lon)&&isvector(data)
                    figure
                    if not(isempty(shape))
                        latmin=min(lat);
                        latmax=max(lat);
                        lonmin=min(lon);
                        lonmax=max(lon);
                        rlat=latmax-latmin;
                        rlon=lonmax-lonmin;
                        mapshow(shape);
                        xlim([lonmin-rlon*0.1,lonmax+rlon*0.1]);
                        ylim([latmin-rlat*0.1,latmax+rlat*0.1]);
                        hold on
                    end
                    minval=nanmin(data);
                    maxval=nanmax(data);
                    for i=1:length(lat)
                        if not(isnan(data(i)))
                            %temp1=(data(i)-minval)/(maxval-minval);
                            %color=[temp1 temp1 temp1];
                            color = stem_misc.get_rainbow_color(data(i),minval,maxval);
                            mapshow(lon(i),lat(i),'DisplayType','point','MarkerFaceColor',color, 'MarkerEdgeColor','k','Marker','o','MarkerSize',5);
                            hold on
                        end
                    end
                else
                    error('lat, lon and data must be either all matrices or all vectors');
                end
            end
        end
        
        function v = triuv(mat)
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
        
        function weights = wendland(DistMat,gamma,resh)
            if gamma<=0
                error('The gamma tapering parameter must be > 0');
            end
            if nargin<3
                resh=0;
            end
            %find the indices of the non-zero elements of DistMat;
            idx=find(DistMat);
            %evaluate the wendland function only on the non-zero elements
            c=DistMat(idx)./gamma;
            %weights=((1-c).^6).*(35*c.^2+18*c+3);
            weights=((1-c).^4).*(4*c+1);
            if resh
                [I,J]=ind2sub(size(DistMat),idx);
                weights=sparse(I,J,weights,size(DistMat,1),size(DistMat,2));
            end
        end

    end
end