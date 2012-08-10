res=st_loo_residual.Y./sqrt(st_loo_residual.kriging_Var_W_bar_hat);
res1=res(1:66,:);
res2=res(67:76,:);
res3=res(77:end,:);
[f1,x1]=ksdensity(res1(:),'width',1,'npoints',1000);
[f2,x2]=ksdensity(res2(:),'width',1,'npoints',1000);
[f3,x3]=ksdensity(res3(:),'width',1,'npoints',1000);
plot(x1,f1,'b');
hold on
plot(x2,f2,'b--');
plot(x3,f3,'b:');
