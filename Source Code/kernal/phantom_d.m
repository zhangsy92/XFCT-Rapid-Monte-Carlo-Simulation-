% Calculate the path length of a particle in a circle region
function dist=phantom_d(xs,ys,xed,yed,xc,yc,rc)
dist=zeros(size(xed));
xsc=xc-xs;
ysc=yc-ys;
xsed=xed-xs;
ysed=yed-ys;
lsc=sqrt(xsc^2+ysc^2);
lsed=sqrt(xsed.^2+ysed.^2);
costl=(xsc*xsed+ysc*ysed)./(lsc.*lsed);
costc=sqrt(lsc^2-rc^2)/lsc;
lef=costl>costc;
dl=sqrt(rc^2-lsc^2.*(1-costl.^2));
d1=lsc*costl-dl;
d2=max(d1,min(lsed,d1+dl*2));
dist(lef)=d2(lef)-d1(lef);
end