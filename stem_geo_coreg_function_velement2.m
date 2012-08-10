%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = stem_geo_coreg_function_velement2(v_element,row,col,v_full,theta,correlation_type,DistMat,var_dims,E_w_y1,sum_Var_w_y1,tapering_par,r)
%only V is estimated

if isempty(tapering_par)
    error('This method is available only with tapering');
end

block_size=5000;

n_var=length(var_dims);

v=v_full;
v(row,col)=v_element;
v(col,row)=v_element;

if min(eig(v))>0
    sigma_W=DistMat;

    for j=1:n_var
        for i=j:n_var
            [B,block_i,block_j] = get_block(var_dims,i,var_dims,j,DistMat);
            sigma_W(block_i,block_j)=v(i,j)*stem_correlation_function(theta,B,correlation_type);
            sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*stem_wendland(DistMat(block_i,block_j),tapering_par);
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
    
    blocks=0:block_size:length(r);
    if blocks(end)<length(r)
        blocks=[blocks,length(r)];
    end
   
    E_w_y1=E_w_y1(r,:);
    sum_Var_w_y1=sum_Var_w_y1(r,r);
    c=chol(sigma_W(r,r));
    f=2*T*sum(log(diag(c)));
    for i=1:length(blocks)-1
        idx=blocks(i)+1:blocks(i+1);
        total=zeros(size(E_w_y1,1),length(idx));
        for t=1:size(E_w_y1,2)
            total=total+E_w_y1(:,t)*E_w_y1(idx,t)';
        end
        total=total+sum_Var_w_y1(:,idx);
        res=sigma_W(r,r)\total;
        res=res((1:length(idx))+blocks(i),:);
        f=f+trace(res);
    end    
else
    f=10^10;
end
%f=T*2*sum(log(diag(chol(sigma_W_core))))+trace(sigma_W_core\U); 