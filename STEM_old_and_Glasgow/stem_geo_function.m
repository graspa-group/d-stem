function f = stem_geo_function(theta,correlation_type,DistMat,U,T)

if strcmp(correlation_type,'exponential')
    sigma_W=exp(-DistMat./theta);
end
if strcmp(correlation_type,'matern')
    sigma_W=stem_matern_function(theta(1),theta(2),DistMat);
    sigma_W(isnan(sigma_W))=1;
end


f=T*sum(log(eig(sigma_W)))+trace(sigma_W\U);
%f=T*2*sum(log(diag(chol(sigma_W))))+trace(sigma_W\U);