%------------------------%
%本程序获取重建所需系统矩阵
%------------------------%
clear
clc
set(0,'defaultfigurecolor','w')

load('uatt_90keV_double.mat')
p=log(max(max(n_detR))./(n_detR));
theta=-(0:1:359)*(pi/180);
DetWidth=12.8;
DistSourceAxis=30;
DistSourceDet=60;
N=70;
PixSize=0.05;
u_map=FANFBP_ED_v2(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize);
p=log(max(max(n_detL))./(n_detL));
u_map2=FANFBP_ED_v2(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize);
u_map=u_map';
u_map2=u_map2';


l=70;
disDec=30;%射线源到模体中心距离
p=0.05;%探测器晶体厚度
h=0.2;%探测器晶体高度
D=0.05;%准直孔宽度
b=4;%摸体中心到孔距离
a=6.35;%探测器中心到孔距离
cij=sparse(360*256*2,l^2);

n=1:l^2;
for t=0:1:359
    t
    u=imrotate(u_map,t,'bilinear');
    u2=imrotate(u_map2,t,'bilinear');
    x1=mod(n-1,l)*0.05-0.05*(l/2-0.5);
    y1=ceil(n/l)*0.05-0.05*(l/2+0.5);
    G=[cos(t/180*pi) sin(t/180*pi);-sin(t/180*pi) cos(t/180*pi)];
    nxy=[x1' y1']*G;
    x=nxy(:,1);
    y=nxy(:,2);
    theta=atan(y./(x+disDec));
    x2=-l/2*0.05;
    y2=y-(x-x2).*tan(theta);
    d=sqrt((x-x2).^2+(y-y2).^2);
    avg_u=zeros(length(x),1);
    
    x3=0;
    y3=-b;
    d3=sqrt((x-x3).^2+(y-y3).^2);
    d23=sqrt((x-x3).^2+(y+y3).^2);
    avg_u3=zeros(length(x),1);
    avg2_u3=zeros(length(x),1);
    
    for m=1:100
        xu=x2+m*(x-x2)/100;xu=ceil(xu./0.05+length(u)/2);
        yu=y2+m*(y-y2)/100;yu=ceil(yu./0.05+length(u)/2);
        xu3=x3+m*(x-x3)/100;xu3=ceil(xu3./0.05+length(u)/2);
        yu3=y3+m*(y-y3)/100;yu3=ceil(yu3./0.05+length(u)/2);
        y2u3=y3+m*(y+y3)/100;y2u3=ceil(y2u3./0.05+length(u)/2);
        for q=1:length(x)
            avg_u(q)=avg_u(q)+u(xu(q),yu(q))/100;
            if yu3(q)>0
                avg_u3(q)=avg_u3(q)+u2(xu3(q),yu3(q))/100;
            end
            if y2u3(q)>0
                avg2_u3(q)=avg2_u3(q)+u2(xu3(q),y2u3(q))/100;
            end
        end
    end
    uatt=exp(-d.*avg_u);
    uatt3=exp(-d3.*avg_u3);
    u2att3=exp(-d23.*avg2_u3);
    
    pinphi1=atan((-D/2-x)./(b+y));
    pinphi2=atan((D/2-x)./(b+y));
    pin2phi1=atan((x-D/2)./(b-y));
    pin2phi2=atan((x+D/2)./(b-y));
    
    k=1:256;k=repmat(k,l^2,1);
    
    detx1=k*p-128.5*p-p/2;
    detx2=k*p-128.5*p+p/2;
    det2x1=-detx1;
    det2x2=-detx2;
    
    x=repmat(x,1,256);
    y=repmat(y,1,256);
    pinphi1=repmat(pinphi1,1,256);
    pinphi2=repmat(pinphi2,1,256);
    pin2phi1=repmat(pin2phi1,1,256);
    pin2phi2=repmat(pin2phi2,1,256);
    
    detphi1=atan((detx1-x)./(a+b+y));
    detphi2=atan((detx2-x)./(a+b+y));
    det2phi1=atan((x-det2x1)./(a+b-y));
    det2phi2=atan((x-det2x2)./(a+b-y));
    
    dphi=max(0,min(detphi2,pinphi2)-max(detphi1,pinphi1));
    d2phi=max(0,min(det2phi2,pin2phi2)-max(det2phi1,pin2phi1));
    
    dh=h./sqrt((a+b+y).^2+(k*p-128.5*p-x).^2);
    d2h=h./sqrt((a+b-y).^2+(-(k*p-128.5*p)-x).^2);
%     cij(t*256+1:t*256+256,:)=dh'.*dphi'.*repmat(uatt',256,1).*repmat(uatt3',256,1);
    
    cij(t*512+1:t*512+256,:)=dh'.*dphi'.*repmat(uatt',256,1).*repmat(uatt3',256,1);
    cij(t*512+257:t*512+512,:)=d2h'.*d2phi'.*repmat(uatt',256,1).*repmat(u2att3',256,1);
    
    
    
end
save projdata\SMsparse_att_Gd_double cij


