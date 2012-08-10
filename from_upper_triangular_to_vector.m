%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vec = from_upper_triangular_to_vector(mat)
 d=size(mat,1);
 vec=zeros(d*(d-1)/2,1);
 
 counter=1;
 for i=1:d-1
     for j=i+1:d
         vec(counter)=mat(i,j);
         counter=counter+1;
     end
 end
 
end