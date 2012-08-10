%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = stem_geo_coreg_function_velement(v_element,row,col,v_full,theta,correlation_type,DistMat,var_dims,U,T,tapering_par,r)
%only V is estimated
n_var=length(var_dims);

v=v_full;
v(row,col)=v_element;
v(col,row)=v_element;

if min(eig(v))>0
    if not(isempty(tapering_par))
        sigma_W=DistMat;
    else
        sigma_W=zeros(sum(var_dims));    
    end
    for j=1:n_var
        for i=j:n_var
            [B,block_i,block_j] = get_block(var_dims,i,var_dims,j,DistMat);
            sigma_W(block_i,block_j)=v(i,j)*stem_correlation_function(theta,B,correlation_type);
            if not(isempty(tapering_par))
                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*stem_wendland(DistMat(block_i,block_j),tapering_par);
            end            
            if (i~=j)
                sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
            end
        end
    end
    if strcmp(correlation_type,'matern')
        for i=1:size(sigma_W,1)
            sigma_W(i,i)=1;
        end
        sigma_W(isnan(sigma_W_core))=0;
    end
    if not(isempty(tapering_par))
        %r = symamd(sigma_W);
        c=chol(sigma_W(r,r));
        f=2*T*sum(log(diag(c)))+trace(chol_solve(c,U(r,r)));
    else
        c=chol(sigma_W);
        f=2*T*sum(log(diag(c)))+trace(chol_solve(c,U));
    end
else
    f=10^10;
end
%f=T*2*sum(log(diag(chol(sigma_W_core))))+trace(sigma_W_core\U); 