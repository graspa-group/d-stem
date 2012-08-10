function res_kalman=get_temporal_residual(res_ols,zsmooth,K)
  res_kalman=zeros(size(res_ols));
  for t=1:length(zsmooth)
     res_kalman(:,t)=res_ols(:,t)-K(:,:,t)*zsmooth(:,t);
  end
end
