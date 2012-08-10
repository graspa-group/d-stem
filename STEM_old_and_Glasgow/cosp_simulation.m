%% pixel vicinato ordine 1
semiside=0.5;
x1=unifrnd(-semiside,semiside,15000000,1);
y1=unifrnd(-semiside,semiside,15000000,1);
x2=unifrnd(-semiside,semiside,15000000,1);
y2=unifrnd(-semiside,semiside,15000000,1);
%0 gradi
c1=[0,0];
c2=[5,5];
%d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
d=sqrt(((c1(1))-(c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
[f,x]=hist(d,100);
hist(d,100);
level1(1).degree=0;
level1(1).x_distance=1;
level1(1).y_distance=0;
level1(1).pixel_side=1;
level1(1).f=f/sum(f);
level1(1).x=x;
%45 gradi
c1=[0,0];
c2=[1,1];
d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
[f,x]=hist(d,100);
level1(2).degree=45;
level1(2).x_distance=1;
level1(2).y_distance=1;
level1(2).pixel_side=1;
level1(2).f=f/sum(f);
level1(2).x=x;

%% pixel vicinato ordine 2
%0 gradi
c1=[0,0];
c2=[0,2];
d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
[f,x]=hist(d,100);
level2(1).degree=0;
level2(1).x_distance=2;
level2(1).y_distance=0;
level2(1).pixel_side=1;
level2(1).f=f/sum(f);
level2(1).x=x;
%30 gradi
c1=[0,0];
c2=[1,2];
d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
[f,x]=hist(d,100);
level2(2).degree=30;
level2(2).x_distance=2;
level2(2).y_distance=1;
level2(2).pixel_side=1;
level2(2).f=f/sum(f);
level2(2).x=x;
%45 gradi
c1=[0,0];
c2=[2,2];
d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
[f,x]=hist(d,100);
level2(3).degree=45;
level2(3).x_distance=2;
level2(3).y_distance=2;
level2(3).pixel_side=1;
level2(3).f=f/sum(f);
level2(3).x=x;

%% pixel vicinato ordine 100
for i=0:100
    i
    %0 gradi
    c1=[0,0];
    c2=[100,i];
    d=sqrt(((x1-c1(1))-(x2-c2(1))).^2+((y1-c1(2))-(y2-c2(2))).^2);
    [f,x]=hist(d,100);
    level3(i+1).degree=rad2deg(atan(c2(2)/c2(1)));
    level3(i+1).x_distance=100;
    level3(i+1).y_distance=i;
    level3(i+1).pixel_side=1;
    level3(i+1).f=f/sum(f);
    level3(i+1).x=x;
end