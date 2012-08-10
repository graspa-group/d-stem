%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = stem_geo_coreg_function_theta2(theta,v_full,correlation_type,DistMat,var_dims,E_w_y1,sum_Var_w_y1,tapering_par,r)
%only theta is estimated
if isempty(tapering_par)
    error('This method is available only with tapering');
end

block_size=2000;

T=size(E_w_y1,2);
n_var=length(var_dims);
v=v_full;

if min(eig(v))>0
    sigma_W=DistMat;

    I=zeros(nnz(DistMat),1);
    J=zeros(nnz(DistMat),1);
    elements=zeros(nnz(DistMat),1);
    idx=0;
    blocks=[0 cumsum(var_dims)];
    for j=1:n_var
        for i=j:n_var
            [B,block_i,block_j] = get_block(var_dims,i,var_dims,j,DistMat);
            corr_result=stem_correlation_function(theta,B,correlation_type);
            weights=stem_wendland(B,tapering_par); %possibile calcolarli una sola volta???
            corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
            l=length(corr_result.I);
            I(idx+1:idx+l)=corr_result.I+blocks(i);
            J(idx+1:idx+l)=corr_result.J+blocks(j);
            elements(idx+1:idx+l)=corr_result.correlation;
            idx=idx+l;
            if not(i==j)
                I(idx+1:idx+l)=corr_result.J+blocks(j);
                J(idx+1:idx+l)=corr_result.I+blocks(i);
                elements(idx+1:idx+l)=corr_result.correlation;
                idx=idx+l;
            end
        end
    end
    sigma_W=sparse(I,J,elements);
    
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