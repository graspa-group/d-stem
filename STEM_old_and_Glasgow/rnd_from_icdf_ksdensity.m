function  y = rnd_from_icdf_ksdensity(f,xi,n)
    u=rand(n,1);
    u=round(u*length(xi))+1;
    %u(u>length(xi))=length(xi);
    y=f(u);
end