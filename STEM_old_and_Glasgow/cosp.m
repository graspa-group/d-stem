clc
clear all
load ../Data/cosp_level_structure.mat
theta=50;
semiside=0.5;

x1=unifrnd(-semiside,semiside,100000,1);
y1=unifrnd(-semiside,semiside,100000,1);
x2=unifrnd(-semiside,semiside,100000,1);
y2=unifrnd(-semiside,semiside,100000,1);
c1.x=0;
c1.y=0;
c2.x=0;
c2.y=1;
d=sqrt(((x1-c1.x)-(x2-c2.x)).^2+((y1-c1.y)-(y2-c2.y)).^2);
%[f,x]=hist(d,100);
%f=f/sum(f);
%corr=exp(-x/theta)*f'
corr=mean(exp(-d/theta))


degree_list=[];
for i=1:length(level3)
    degree_list(i)=level3(i).degree;
end

tic
angle=rad2deg(atan(abs(c1.x-c2.x)/abs(c2.y-c2.y)));
if angle>45
    angle=90-angle;
end
delta_angle=abs(degree_list-angle);
[value,index]=min(delta_angle);
p1a.x=c1.x-semiside;
p1a.y=c1.y-semiside;
p1b.x=c1.x+semiside;
p1b.y=c1.y+semiside;
p2a.x=c2.x+semiside;
p2a.y=c2.y+semiside;
p2b.x=c2.x-semiside;
p2b.y=c2.y-semiside;

d1=(p1a.x-p2a.x)^2+(p1a.y-p2a.y)^2;
d2=(p1a.x-p2b.x)^2+(p1a.y-p2b.y)^2;
d3=(p1b.x-p2b.x)^2+(p1b.y-p2b.y)^2;
dmin=sqrt(min([d1,d2,d3]));
dmax=sqrt(max([d1,d2,d3]));

if max(abs(c1.x-c2.x),abs(c1.y-c2.y))==1
    %level1

else
    if max(abs(c1.x-c2.x),abs(c1.y-c2.y))==2
        %level2
    else
        %level3
        x=dmin:(dmax-dmin)/(length(level3(1).x)-1):dmax;
        corr_approx=exp(-x/theta)*level3(index).f';
    end
end
toc
abs(corr-corr_approx)