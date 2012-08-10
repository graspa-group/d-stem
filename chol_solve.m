%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        [I,J]=find(abs(y)>1e-5);
        L=sub2ind(size(y),I,J);
        elements=y(L);
        y=sparse(I,J,elements,size(y,1),size(y,2));
        res=c\y;
        [I,J]=find(abs(res)>1e-5);
        L=sub2ind(size(res),I,J);
        elements=res(L);
        res=sparse(I,J,elements,size(res,1),size(res,2));     
    end

end