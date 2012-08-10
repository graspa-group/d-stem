%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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