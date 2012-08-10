%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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