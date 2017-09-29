%-----------------------------%
%Rapid Monte Carlo simulation of Au x-ray fluorescent photon detection in XFCT system
%X-ray Source: I=100mA V=150kVp  filter:2mm Cu
%Materials: Au H2O PMMA
%Acquisition time per angle: 1s
%2017-04-01
%-----------------------------%
clear;
clc;
load('data\spec_150kVp_2mmCu.mat');
load('data\Um.mat');

%Initialization of fluorescent photon count
n_det1=zeros(256,360);
n_det2=zeros(256,360);

%% Parameters
current=100;                   %Current
t=1;                           %Acquisition time per angle
p=0.05;                        %width of the detector pixel
keel=0.02;                     %keel-edge
h=0.2;                         %height of the pinhole
D=0.05;                        %width of the pinhole
b=4;                           %Distance between rotation center and pinhole
a=6.35;                        %Distance between detector center and pinhole
disDec=30;                     %Distance between x-ray source and rotation center
n=5e6;                         %Sample number
R=1.75;                        %Radius of phantom
r=0.26;                        %Radius of solution

rou_Au=[0.001;0.002;0.003;0.004];  %Concentration of solutions (w/v)

% Cross section of xrf photons
uk1li=exp(interp1(eli,log(uli),68.794));
uk1pmma=exp(interp1(epmma,log(upmma),68.794))*1.18;
uk1au=exp(interp1(eau,log(uau),68.794));
uk2li=exp(interp1(eli,log(uli),66.98));
uk2pmma=exp(interp1(epmma,log(upmma),66.98))*1.18;
uk2au=exp(interp1(eau,log(uau),66.98));

%% Particle sampling
phi=asin(R/disDec);
Ep=80+rand(n,1)*70;             %Energy sampling
Np=interp1(E,N,Ep);             %Spectrum interpolation
thetap=-phi+rand(n,1)*2*phi;    %Incident angle sampling

%Coordinate of insertions
x_Au=[1;0.5;-0.5;-1];
y_Au=[0;0.5*sqrt(3);0.5*sqrt(3);0];

%% Fluorescenct photon count
nk=0;
for k=0:1/180*pi:2*pi-1/180*pi
    nk=nk+1; %Projection number
    fprintf('projection %d\n',nk)
    
    Ixf=[];X_xf=[];Y_xf=[];
    
    %Coordinate of insertions after rotation
    G=[cos(k) sin(k);-sin(k) cos(k)];
    nxy=[x_Au y_Au]*G;
    X_Au=nxy(:,1);
    Y_Au=nxy(:,2);
    
    for n_Au=1:4 %Calculate xrf photons in different insertions
        tt=atan(Y_Au(n_Au)/(X_Au(n_Au)+disDec));dt1=asin(r/sqrt(Y_Au(n_Au)^2+(X_Au(n_Au)+disDec)^2));
        isxf=abs(thetap-tt)<dt1;
        I1=Np(isxf);E1=Ep(isxf);C1=thetap(isxf);
        dm1=2*sqrt(r^2-(sqrt(Y_Au(n_Au)^2+(X_Au(n_Au)+disDec)^2).*sin(C1-tt)).^2); %Max path length of particles in the insertion
        
        %Total cross section of the insertion
        Miu1=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(n_Au))*(1+rou_Au(n_Au));
        
        L1=-log(rand(length(I1),1))./Miu1;%Sample path length of particles before interaction
        isxf=L1<dm1;  %Ignore path length sampled which is longer than maximum path length in the insertion
        
        I1=I1(isxf);L1=L1(isxf);C1=C1(isxf);dm1=dm1(isxf);E1=E1(isxf);Miu1=Miu1(isxf);
        
        %Calculate the coordinate where xrf photons emitted from
        dl1=sqrt(Y_Au(n_Au)^2+(X_Au(n_Au)+disDec)^2).*cos(C1-tt)-disDec*cos(C1);
        d1=sqrt(R^2-(disDec*sin(C1)).^2)+dl1-dm1/2;
        dist=disDec*cos(C1)+dl1-dm1/2+L1;
        y_xf=dist.*sin(C1);
        x_xf=dist.*cos(C1)-disDec;
        
        %Calculate the path length and attenuation coefficient in the phantom
        l_PMMA=phantom_d(-disDec,0,x_xf,y_xf,0,0,R);
        l_Au=zeros(length(x_xf),4);
        miu_Au=zeros(length(x_xf),4);
        for n_at=1:4
            l_Au(:,n_at)=phantom_d(-disDec,0,x_xf,y_xf,X_Au(n_at),Y_Au(n_at),r);
            miu_Au(:,n_at)=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(n_at))*(1+rou_Au(n_at));
        end
        
        %Calculate the number of xrf photons
        I1=I1.*exp(-exp(interp1(epmma,log(upmma),E1))*1.18.*(l_PMMA-sum(l_Au,2))-sum(l_Au.*miu_Au,2)+l_Au(n_Au).*miu_Au(n_Au))...
            .*exp(interp1(eau,log(auphk),E1))*rou_Au(n_Au)./Miu1;
        
        Ixf=[Ixf;I1];
        X_xf=[X_xf;x_xf];
        Y_xf=[Y_xf;y_xf];
    end

    %Detector 1
    x=X_xf;
    y=Y_xf;
    %Equivalent height of beam
    mh=h*(0.5+0.5*(b+y-keel/2)./a)./(0.5+(b+y-keel/2)./a);
    %Emission angle sampling
    phi1=atan((D/2-x)./(b+y-keel/2));phi2=atan((-x-D/2)./(b+y-keel/2));
    phiflu=phi1+rand(length(x),1).*(phi2-phi1);
    %Pixel number
    numdet=floor(((a+b+y).*tan(phiflu)+x)*20+129);
    %Solid angle
    sangle=(b+y-keel/2)./(a+b+y).*mh.*D.*cos(phiflu)./(4*pi*((b+y-keel/2)./cos(phiflu)).^2);
    Id=Ixf.*sangle;
    %Attenuation of xrf
    lxf_PMMA=phantom_d(0,-b,x,y,0,0,R);
    lxf_Au=zeros(length(x),4);
    for n_at=1:4
        lxf_Au(:,n_at)=phantom_d(0,-b,x,y,X_Au(n_at),Y_Au(n_at),r);
    end
    I1=Id*0.482.*exp(-sum(lxf_Au*diag((uk1li+uk1au*rou_Au).*(1+rou_Au)),2)-uk1pmma*1.18.*(lxf_PMMA-sum(lxf_Au,2)));
    I2=Id*0.267.*exp(-sum(lxf_Au*diag((uk2li+uk2au*rou_Au).*(1+rou_Au)),2)-uk2pmma*1.18.*(lxf_PMMA-sum(lxf_Au,2)));
    %Count
    for i=1:length(Ixf)
        n_det1(numdet(i),nk)=n_det1(numdet(i),nk)+I1(i)+I2(i);
    end
    
    %Detector 2
    x=-X_xf;
    y=-Y_xf;
    %Equivalent height of beam
    mh=h*(0.5+0.5*(b+y-keel/2)./a)./(0.5+(b+y-keel/2)./a);
    %Emission angle
    phi1=atan((D/2-x)./(b+y-keel/2));phi2=atan((-x-D/2)./(b+y-keel/2));
    phiflu=phi1+rand(length(x),1).*(phi2-phi1);
    %Pixel number
    numdet=floor(((a+b+y).*tan(phiflu)+x)*20+129);
    %Solid angle
    sangle=(b+y-keel/2)./(a+b+y).*mh.*D.*cos(phiflu)./(4*pi*((b+y-keel/2)./cos(phiflu)).^2);
    Id2=Ixf.*sangle;
    %Attenuation of xrf
    lxf_PMMA=phantom_d(0,-b,x,y,0,0,R);
    lxf_Au=zeros(length(x),4);
    for n_at=1:4
        lxf_Au(:,n_at)=phantom_d(0,-b,x,y,X_Au(n_at),Y_Au(n_at),r);
    end
    I1=Id*0.482.*exp(-sum(lxf_Au*diag((uk1li+uk1au*rou_Au).*(1+rou_Au)),2)-uk1pmma*1.18.*(lxf_PMMA-sum(lxf_Au,2)));
    I2=Id*0.267.*exp(-sum(lxf_Au*diag((uk2li+uk2au*rou_Au).*(1+rou_Au)),2)-uk2pmma*1.18.*(lxf_PMMA-sum(lxf_Au,2)));
    %Count
    for i=1:length(Ixf)
        n_det2(numdet(i),nk)=n_det2(numdet(i),nk)+I1(i)+I2(i);
    end
end

%% Total dose calculation
area=2*disDec*tan(phi)*p*(0.5*a+b)/(0.5*a)/(0.2^2);%dose area
N=N.*area*current*t;% total number of incident photons
%Absolute dose conversion
sig1=sum(N(find(E==80):find(E==150)))/sum(Np);
n_det1=n_det1*sig1;
n_det2=n_det2*sig1;
n_det1=random('poisson',n_det1);
n_det2=random('poisson',n_det2);
%% Save
save projdata\xrf_Au_150keV_double n_det1 n_det2