%-----------------------------%
%Calculation of attenuation map for Gd
%X-ray Source: I=100mA V=90kVp  filter:0.3mm Cu
%Materials: Gd H2O PMMA
%Acquisition time per angle: 1s
%2017-04-01
%-----------------------------%clear;
clc;
load('data\spec_150kVp_2mmCu.mat');
load('data\Um.mat');

%Initialization of CT count
n_detL=zeros(256,360);
n_detR=zeros(256,360);

%% Parameters
current=100;                   %Current
t=1;                           %Acquisition time per angle
p=0.05;                        %width of the detector pixel
h=0.2;                         %height of the pinhole
disDec=30;                     %Distance between x-ray source and rotation center
disSD=60;                      %Distance between x-ray source and detector
R=1.75;                        %Radius of phantom
r=0.26;                        %Radius of solution
rou_Au=[0.001;0.002;0.003;0.004];  %Concentration of solutions (w/v)

%% Particle sampling
tphi=-128*p+p/2:p:128*p-p/2;
dphi=atan(tphi/disSD);
lphi=atan((tphi-p/2)/disSD);
rphi=atan((tphi+p/2)/disSD);
Dphi=atan((tphi+p/2)/disSD)-atan((tphi-p/2)/disSD);
dEpL=65:0.1:70;
dEpR=85:0.1:90;
thetap=zeros(length(dphi)*length(dEpL),1);
EpL=thetap;NpL=thetap;
EpR=thetap;NpR=thetap;
n=length(dphi)*length(dEpL);
for ii=1:length(dphi)*length(dEpL)
    thetap(ii)=lphi(ceil(ii/length(dEpL)))+rand(1)*(rphi(ceil(ii/length(dEpL)))-lphi(ceil(ii/length(dEpL))));
    EpL(ii)=dEpL(mod(ii-1,length(dEpL))+1);
    EpR(ii)=dEpR(mod(ii-1,length(dEpR))+1);
    NpL(ii)=interp1(E,N,EpL(ii))*Dphi(ceil(ii/length(dEpL)));
    NpR(ii)=interp1(E,N,EpR(ii))*Dphi(ceil(ii/length(dEpR)));
end

%Coordinate of insertions
x_Au=[1;0.5;-0.5;-1];
y_Au=[0;0.5*sqrt(3);0.5*sqrt(3);0];

x_det=(disSD-disDec)*ones(length(thetap),1);
y_det=disSD*tan(thetap);

%% CT data calculation
%air data
n_airL=zeros(256,360);
n_airR=zeros(256,360);
    for i=1:n
        j=floor(disSD/p*tan(thetap(i))+129);
        n_airL(j,:)=n_airL(j,:)+NpL(i);
        n_airR(j,:)=n_airR(j,:)+NpR(i);
    end

%projection data
nk=0;
for k=0:1/180*pi:2*pi-1/180*pi
    nk=nk+1;
    fprintf('projection %d\n',nk)
    %Coordinate of insertions after rotation
    G=[cos(k) sin(k);-sin(k) cos(k)];
    nxy=[x_Au y_Au]*G;
    X_Au=nxy(:,1);
    Y_Au=nxy(:,2);
    
    %Calculate the path length and attenuation coefficient in the phantom
    l_PMMA=phantom_d(-disDec,0,x_det,y_det,0,0,R);
    miu_PMMAL=1.18*exp(interp1(epmma,log(upmma),EpL));
    miu_PMMAR=1.18*exp(interp1(epmma,log(upmma),EpR));
    
    l_Au=zeros(length(x_det),4);
    miu_AuL=zeros(length(x_det),4);
    miu_AuR=zeros(length(x_det),4);
    for n_at=1:4
        l_Au(:,n_at)=phantom_d(-disDec,0,x_det,y_det,X_Au(n_at),Y_Au(n_at),r);
        miu_AuL(:,n_at)=(exp(interp1(eli,log(uli),EpL))+exp(interp1(eau,log(uau),EpL))*rou_Au(n_at))*(1+rou_Au(n_at));
        miu_AuR(:,n_at)=(exp(interp1(eli,log(uli),EpR))+exp(interp1(eau,log(uau),EpR))*rou_Au(n_at))*(1+rou_Au(n_at));
    end
    
    %Calculate count after attenuation
    IL=NpL.*exp(-sum(l_Au.*miu_AuL,2)-(l_PMMA-sum(l_Au,2)).*miu_PMMAL);
    IR=NpR.*exp(-sum(l_Au.*miu_AuR,2)-(l_PMMA-sum(l_Au,2)).*miu_PMMAR);
    
    %Pixel number
    for i=1:n
        j=floor(disSD/p*tan(thetap(i))+129);
        n_detL(j,nk)=n_detL(j,nk)+IL(i);
        n_detR(j,nk)=n_detR(j,nk)+IR(i);
    end

end

%% Total dose calculation
area=256*p*h*(100/disSD)^2;%dose area
N=N.*area*current*t;% total number of incident photons
%Absolute dose conversion
sig1=sum(N(find(E==65):find(E==70)))/sum(NpL);
sig2=sum(N(find(E==85):find(E==90)))/sum(NpR);
n_detL=n_detL*sig1;
n_detR=n_detR*sig2;
n_detL=random('poisson',n_detL);
n_detR=random('poisson',n_detR);
n_airL=n_airL*sig1;
n_airR=n_airR*sig2;
n_airL=random('poisson',n_airL);
n_airR=random('poisson',n_airR);

pr=-log(n_detR./n_airR);
pl=-log(n_detL./n_airL);
theta=-(0:1:359)*(pi/180);
DetWidth=12.8;
N=70;
PixSize=0.05;
% Attenuation map reconstruction using filtered back projection algorithm
u_map1=FANFBP_ED_v2(pr,theta,DetWidth,disDec,disSD,N,PixSize)';
u_map2=FANFBP_ED_v2(pl,theta,DetWidth,disDec,disSD,N,PixSize)';

%% Save
save projdata\uatt_150keV_double u_map1 u_map2

