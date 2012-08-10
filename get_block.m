%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [B,block_i,block_j] = get_block(dim_r,i,dim_c,j,A)
% [B,block_i,block,j] = get_block(dim_r,i,dim_c,j,A)
%
% rende B=A(block_i,block,j)
% dim_r = struttura a blocchi delle righe
% dim_c = struttura a blocchi delle colonne

rr=cumsum(dim_r);
block_i = [rr(i)-dim_r(i)+1:rr(i)];
if nargin<=2
    % rende solo block_i
    B=block_i;
    block_i=[];
    block_j=[];
    return
end

cc=cumsum(dim_c);
block_j = [cc(j)-dim_c(j)+1:cc(j)];
if nargin<=4
    B=block_i;
    block_i=block_j;
    block_j=[];
    return
end

B=A(block_i,block_j);


end