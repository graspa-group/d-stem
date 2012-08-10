%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dist = get_distance(lat1,lon1,lat2,lon2)
 dlat=(lat1-lat2)/360*2*pi;
 dlong=(lon1-lon2)/360*2*pi;
 lat1=lat1/360*2*pi;
 lat2=lat2/360*2*pi;
 %lon1=lon1/360*2*pi;
 %lon2=lon2/360*2*pi;
 a=sin(dlat/2).^2+cos(lat1).*cos(lat2).*sin(dlong/2).^2;
 c=2*atan2(a.^0.5,(1-a).^0.5);
 dist=6371*c;
end