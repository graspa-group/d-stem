%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weights = stem_wendland(DistMat,gamma,resh)
    if gamma<=0
        error('The gamma tapering parameter must be > 0');
    end
    if nargin<3
        resh=0;
    end
    %find the indices of the non-zero elements of DistMat;
    idx=find(DistMat);
    %evaluate the wendland function only on the non-zero elements
    c=DistMat(idx)./gamma;
    %weights=((1-c).^6).*(35*c.^2+18*c+3);
    weights=((1-c).^4).*(4*c+1);
    if resh
        [I,J]=ind2sub(size(DistMat),idx);
        weights=sparse(I,J,weights,size(DistMat,1),size(DistMat,2));
    end
end