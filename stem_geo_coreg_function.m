function f = stem_geo_coreg_function(theta_and_v,correlation_type,DistMat,var_dims,U,T,cosp_C,cosp_theta_vector)

n_var=length(var_dims);

if strcmp(correlation_type,'exponential')
    theta_coreg=theta_and_v(1);
    v=from_vector_to_symmetric_matrix(theta_and_v(2:end));
end
if strcmp(correlation_type,'matern')
    theta_coreg=theta_and_v(1:2);
    v=from_vector_to_symmetric_matrix(theta_and_v(3:end));
end

if min(eig(v))>0
    sigma_W_core=zeros(sum(var_dims));    
    for j=1:n_var
        for i=j:n_var
            [B,block_i,block_j] = get_block(var_dims,i,var_dims,j,DistMat);
            if strcmp(correlation_type,'exponential')
                sigma_W_core(block_i,block_j)=v(i,j)*exp(-B./theta_coreg);
            end
            if strcmp(correlation_type,'matern')
                if not(v(i,j)==0)
                    sigma_W_core(block_i,block_j)=v(i,j)*stem_matern_function(theta_coreg(1),theta_coreg(2),B);
                else
                    sigma_W_core(block_i,block_j)=zeros(length(block_i),length(block_j));
                end
            end
            if (i~=j)
                sigma_W_core(block_j,block_i)=sigma_W_core(block_i,block_j)';
            end
        end
    end
    if strcmp(correlation_type,'matern')
        for i=1:size(sigma_W_core,1)
            sigma_W_core(i,i)=1;
        end
    end
    sigma_W_core(isnan(sigma_W_core))=0;
    if isempty(cosp_C)==0
        if (theta_coreg<max(cosp_theta_vector))&&(theta_coreg>min(cosp_theta_vector))
            for l=1:length(cosp_theta_vector)-1
                if (theta_coreg>=cosp_theta_vector(l))&&(theta_coreg<cosp_theta_vector(l+1))
                    idx=l;
                end
            end
            r=(theta_coreg-cosp_theta_vector(idx))/(cosp_theta_vector(idx+1)-cosp_theta_vector(idx));
            correction=cosp_C(:,:,idx)+(cosp_C(:,:,idx+1)-cosp_C(:,:,idx))*r;
            sigma_W_core=sigma_W_core.*correction;
            if min(eig(sigma_W_core))<0
                error('Matrix non SDP after cosp correction in numerical optimization!');
            end
        end
    end
    f=T*sum(log(eig(sigma_W_core)))+trace(sigma_W_core\U);
else
    f=10^10;
end
%f=T*2*sum(log(diag(chol(sigma_W_core))))+trace(sigma_W_core\U); 