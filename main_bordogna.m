%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all

pathparallel='/opt/matNfs/';

flag_parallel=0;
flag_tapered=1;

if flag_tapered
    load ../Data/st_model_small_tapered_time
else
    load ../Data/st_model_small_nottapered
end

if flag_parallel
    st_model.EM_estimate(0.001,200,'single',pathparallel);
else
    st_model.EM_estimate(0.001,30,'single');
end
st_model.set_Hessian();