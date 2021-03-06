% Calculate the maximum path length of a particle in a circle region
function dist=c_csline(xs,ys,xed,yed,xc,yc,rc)
l=length(xed);
dist=zeros(l,1);
dl=zeros(l,1);
xsc=xc-xs;
ysc=yc-ys;
xsed=xed-xs;
ysed=yed-ys;
lsc=sqrt(xsc^2+ysc^2);
lsed=sqrt(xsed.^2+ysed.^2);
costl=(xsc*xsed+ysc*ysed)./(lsc.*lsed);
costc=sqrt(lsc^2-rc^2)/lsc;
lef=costl>costc;
dl(lef)=sqrt(rc^2-lsc^2.*(1-costl(lef).^2));
d1=lsc*costl(lef)-dl(lef);
d2=max(d1,min(lsed(lef),d1+dl(lef)*2));
dist(lef)=d2-d1;
dist(dist>0 & dist<2*dl)=2*dl(dist>0 & dist<2*dl);
end