%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corr = stem_correlation_function(theta,DistMat,type)
if strcmp(type,'exponential')
    %theta(1)=theta
    if not(isscalar(theta))
        error('Theta must be a scalar');
    end
    if not(issparse(DistMat))
        corr=exp(-DistMat/theta);
    else
        %note that the diagonal of the variance-covariance matrices and
        %possibly some elements of the cross-covariance matrices are equal
        %to eps instead of zero. However, exp(eps)=1.
        idx=find(DistMat);
        correlation=exp(-DistMat(idx)/theta);
        [I,J]=ind2sub(size(DistMat),idx);
        corr.I=I;
        corr.J=J;
        corr.correlation=correlation;
    end
else
    %thata(1)=alpha
    %theta(2)=nu
    if issparse(DistMat)
        error('Sparse matrices not yet supported with Matern correlation function');
    end
    if not(size(theta,2)==2)
        error('Theta must be a 2x1 vector');
    end
    corr=1/(2^(theta(2)-1)*gamma(theta(2)))*((DistMat/theta(1)).^theta(2)).*besselk(theta(2),DistMat/theta(1));
end
end