%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weights = stem_wendland(DistMat,gamma)
    if gamma<=0
        error('The gamma tapering parameter must be > 0');
    end
    %find the indices of the non-zero elements of DistMat;
    idx=find(DistMat);
    %evaluate the wendland function only on the non-zero elements
    c=DistMat(idx)./gamma;
    %weights=((1-c).^6).*(35*c.^2+18*c+3);
    weights=((1-c).^4).*(4*c+1);
end