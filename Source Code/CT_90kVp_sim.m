%-----------------------------%
%Calculation of attenuation map for Gd
%X-ray Source: I=100mA V=90kVp  filter:0.3mm Cu
%Materials: Gd H2O PMMA
%Acquisition time per angle: 1s
%2017-04-01
%-----------------------------%
clear;clc;load('spec_90kVp_0.3mmCu.mat');load('Um.mat');
I=100;t=1;
p=0.05;%探测器晶体厚度
h=0.2;%探测器晶体高度
disDec=30;%射线源到模体中心距离
disSD=60;%射线源到探测器距离
R=1.75;%模体半径
r=0.26;%金溶液半径
rou=0;
rou1=0.001;%金溶液质量百分比
rou2=0.002;
rou3=0.003;
rou4=0.004;

phi=atan(6.4/disSD);
area=256*p*h*(100/disSD)^2;%等效射线面积
N=N.*area*I*t;%总粒子数

tphi=-128*p+p/2:p:128*p-p/2;
dphi=atan(tphi/disSD);
lphi=atan((tphi-p/2)/disSD);
rphi=atan((tphi+p/2)/disSD);
Dphi=atan((tphi+p/2)/disSD)-atan((tphi-p/2)/disSD);
dEpL=40:0.01:45;
dEpR=55:0.01:60;
thetap=zeros(length(dphi)*length(dEpL),1);
EpL=thetap;NpL=thetap;
EpR=thetap;NpR=thetap;
n=length(dphi)*length(dEpL);
for ii=1:length(dphi)*length(dEpL)
%     thetap(ii)=dphi(ceil(ii/length(dEpL)));
    thetap(ii)=lphi(ceil(ii/length(dEpL)))+rand(1)*(rphi(ceil(ii/length(dEpL)))-lphi(ceil(ii/length(dEpL))));
    EpL(ii)=dEpL(mod(ii-1,length(dEpL))+1);
    EpR(ii)=dEpR(mod(ii-1,length(dEpR))+1);
    NpL(ii)=interp1(E,N,EpL(ii))*Dphi(ceil(ii/length(dEpL)));
    NpR(ii)=interp1(E,N,EpR(ii))*Dphi(ceil(ii/length(dEpR)));
end
%模拟粒子数补正
sig1=sum(N(32:37))/sum(NpL);
sig2=sum(N(47:52))/sum(NpR);
sum(N)
n_detL=zeros(256,360);
n_detR=zeros(256,360);


x1=1;y1=0;x2=0.5;y2=0.5*sqrt(3);x3=-0.5;y3=0.5*sqrt(3);x4=-1;y4=0;%金溶液圆柱圆心位置

dcen=abs(disDec*sin(thetap));
dto=2*(R^2-dcen.^2).^0.5;
nk=0;
for k=0:1/180*pi:2*pi-1/180*pi
    nk=nk+1
    G=[cos(k) sin(k);-sin(k) cos(k)];
    nxy=[x1 y1]*G;X1=nxy(1);Y1=nxy(2);
    nxy=[x2 y2]*G;X2=nxy(1);Y2=nxy(2);
    nxy=[x3 y3]*G;X3=nxy(1);Y3=nxy(2);
    nxy=[x4 y4]*G;X4=nxy(1);Y4=nxy(2);%旋转后金溶液圆柱圆心位置
    tt1=atan(Y1/(X1+disDec));dt1=asin(r/sqrt(Y1^2+(X1+disDec)^2));
    tt2=atan(Y2/(X2+disDec));dt2=asin(r/sqrt(Y2^2+(X2+disDec)^2));
    tt3=atan(Y3/(X3+disDec));dt3=asin(r/sqrt(Y3^2+(X3+disDec)^2));
    tt4=atan(Y4/(X4+disDec));dt4=asin(r/sqrt(Y4^2+(X4+disDec)^2));%金溶液圆柱圆心与射线夹角
    
    dcen=abs(disDec*sin(thetap));
    %判断是否穿过模体
    dto=2*max(R^2-dcen.^2,0).^0.5;
    dm1=2*sqrt(max(0,r^2-(sqrt(Y1^2+(X1+disDec)^2).*sin(thetap-tt1)).^2));
    dm2=2*sqrt(max(0,r^2-(sqrt(Y2^2+(X2+disDec)^2).*sin(thetap-tt2)).^2));
    dm3=2*sqrt(max(0,r^2-(sqrt(Y3^2+(X3+disDec)^2).*sin(thetap-tt3)).^2));
    dm4=2*sqrt(max(0,r^2-(sqrt(Y4^2+(X4+disDec)^2).*sin(thetap-tt4)).^2));
    %计算衰减
    IL=NpL.*exp(-dm4.*(exp(interp1(eli,log(uli),EpL))+exp(interp1(eau,log(uau),EpL))*rou4)*(1+rou4)...
        -dm3.*(exp(interp1(eli,log(uli),EpL))+exp(interp1(eau,log(uau),EpL))*rou3)*(1+rou3)...
        -dm2.*(exp(interp1(eli,log(uli),EpL))+exp(interp1(eau,log(uau),EpL))*rou2)*(1+rou2)...
        -dm1.*(exp(interp1(eli,log(uli),EpL))+exp(interp1(eau,log(uau),EpL))*rou1)*(1+rou1)...
        +(dm1+dm2+dm3+dm4-dto).*1.18.*exp(interp1(epmma,log(upmma),EpL)));
    IR=NpR.*exp(-dm4.*(exp(interp1(eli,log(uli),EpR))+exp(interp1(eau,log(uau),EpR))*rou4)*(1+rou4)...
        -dm3.*(exp(interp1(eli,log(uli),EpR))+exp(interp1(eau,log(uau),EpR))*rou3)*(1+rou3)...
        -dm2.*(exp(interp1(eli,log(uli),EpR))+exp(interp1(eau,log(uau),EpR))*rou2)*(1+rou2)...
        -dm1.*(exp(interp1(eli,log(uli),EpR))+exp(interp1(eau,log(uau),EpR))*rou1)*(1+rou1)...
        +(dm1+dm2+dm3+dm4-dto).*1.18.*exp(interp1(epmma,log(upmma),EpR)));
    for i=1:n
        j=floor(1200*tan(thetap(i))+129);%光子落在第j个探测器上
        n_detL(j,nk)=n_detL(j,nk)+IL(i);
        n_detR(j,nk)=n_detR(j,nk)+IR(i);
    end
    n_detL(:,nk)=n_detL(:,nk)*sig1;
    n_detR(:,nk)=n_detR(:,nk)*sig2;
end

p=log((n_detL)./max(max(n_detL)))-log((n_detR)./max(max(n_detR)));
imtool(p,[])
save uatt_90keV_double n_detL n_detR

