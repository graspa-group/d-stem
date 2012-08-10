%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mat = from_vector_to_symmetric_matrix(vec)
 d=(1+sqrt(1+8*length(vec)))/2;
 mat=eye(d);
 
 counter=1;
 for i=1:d-1
     for j=i+1:d
         mat(i,j)=vec(counter);
         mat(j,i)=mat(i,j);
         counter=counter+1;
     end
 end
 
end