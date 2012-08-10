%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = stem_time(s)
    h=floor(s/3600);
    s=s-h*3600;
    m=floor(s/60);
    s=s-m*60;
    if h>0
        str_h=[num2str(h),'h '];
    else
        str_h=[];
    end
    if m>0
        str_m=[num2str(m),'m '];
    else
        str_m=[];
    end
    str_s=[num2str(s,'%5.2f'),'s'];
    str=[str_h,str_m,str_s];
end