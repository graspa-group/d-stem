%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Z = ar1_sim(G,sigma_eta,T,mu0,sigma_0)

Z=zeros(length(G),T+1);
Z0=mvnrnd(mu0,sigma_0)';
Z(:,1)=Z0;
for t=2:T+1
    Z(:,t)=G*Z(:,t-1)+mvnrnd(zeros(length(G),1),sigma_eta)';
end
Z=Z(:,2:T+1);

end



