function absfd = abs(fd)
coef=getcoef(fd);
coef=abs(coef);
absfd =putcoef(fd,coef);
end