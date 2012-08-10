function corr = stem_matern_function(alpha,nu,DistMat)
 corr=1/(2^(nu-1)*gamma(nu))*((DistMat/alpha).^nu).*besselk(nu,DistMat/alpha);
end