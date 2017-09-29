%------------------------%
%System matrix calculation
%------------------------%
clear
clc

load('projdata\uatt_150keV_double.mat')

%% Parameters
l=70;
disDec=30;%Distance between x-ray source and rotation center
p=0.05;%width of the detector pixel
h=0.2;%height of the pinhole
D=0.05;%width of the pinhole
b=4;%Distance between rotation center and pinhole
a=6.35;%Distance between detector center and pinhole
cij=sparse(360*256*2,l^2);

%% System matrix calculation
n=1:l^2;
for t=0:1:359
    fprintf('projection %d\n',t+1)
    
    %Rotate attenuation map to fit current projection angle
    u=imrotate(u_map1,t,'bilinear');
    u2=imrotate(u_map2,t,'bilinear');
    
    % Coordinate of starting point and end point of particles
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
    
    %Attenuation in the phantom
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
    
    uatt=exp(-d.*avg_u);         %Attenuation of incident photons
    uatt3=exp(-d3.*avg_u3);      %Attenuation of fluorescent photons
    u2att3=exp(-d23.*avg2_u3);   %Attenuation of fluorescent photons
    
    %Solid angle from emitted point to detector pixel
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
    
    cij(t*512+1:t*512+256,:)=dh'.*dphi'.*repmat(uatt',256,1).*repmat(uatt3',256,1);
    cij(t*512+257:t*512+512,:)=d2h'.*d2phi'.*repmat(uatt',256,1).*repmat(u2att3',256,1);
end
save projdata\SMsparse_att_Au_double cij


