function is_ok = isnotnan(x)
%
% is_ok = isnotnan(x)
%
% rende "true" se x è not.nan
%
% by AF'04 - 08

k=size(x,2);
for j=1:k
is_ok(:,j) = (~isnan(x(:,j)));
end