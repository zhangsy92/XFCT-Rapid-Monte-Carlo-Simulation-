%-----------------------------%
%Rapid Monte Carlo simulation of Compton scattering detection in XFCT system
%X-ray Source: I=100mA V=150kVp  filter:2mm Cu
%Materials: Au H2O PMMA
%Acquisition time per angle: 1s
%2017-04-01
%-----------------------------%
clear;clc;
load('data\spec_150kVp_2mmCu.mat');
load('data\Um.mat');

%Initialization of scattered photon count
nsca1=zeros(3,256,360);
nsca2=zeros(3,256,360);
%% Parameters
current=100;                         %Current
t=1;                           %Acquisition time per angle
p=0.05;                        %width of the detector pixel
keel=0.02;                     %keel-edge
h=0.2;                         %height of the pinhole
D=0.05;                        %width of the pinhole
b=4;                           %Distance between rotation center and pinhole
a=6.35;                        %Distance between detector center and pinhole
disDec=30;                     %Distance between x-ray source and rotation center
n=5e6;                         %Sample capacity
R=1.75;                        %Radius of phantom
r=0.26;                        %Radius of solution
rou_Au=[0.001;0.002;0.003;0.004];  %Concentration of solution (w/v)

%% Particle sampling
phi=asin(R/disDec);
Ep=60+rand(n,1)*90;           %Energy sampling
Np=interp1(E,N,Ep);           %Spectrum interpolation
thetap=-phi+rand(n,1)*2*phi;  %Incident angle sampling

%Coordinate of insertions
x_Au=[1;0.5;-0.5;-1];
y_Au=[0;0.5*sqrt(3);0.5*sqrt(3);0];

%% Scattered photon count
miu=exp(interp1(epmma,log(upmma),Ep))*1.18; %Attenuation coefficient in the phantom
nk=0;
for k=0:1/180*pi:2*pi-1/180*pi
    nk=nk+1; %Projection number
    fprintf('projection %d\n',nk)
    
    Isca=[];C=[];DL=[];DM=[];LL=[];Es=[];
    
    %Coordinate of insertions after rotation
    G=[cos(k) sin(k);-sin(k) cos(k)];
    nxy=[x_Au y_Au]*G;
    X_Au=nxy(:,1);
    Y_Au=nxy(:,2);
    
    for n_Au=1:4 %Calculate scattered photons in different insertions
        tt=atan(Y_Au(n_Au)/(X_Au(n_Au)+disDec));dt1=asin(r/sqrt(Y_Au(n_Au)^2+(X_Au(n_Au)+disDec)^2));
        isxf=abs(thetap-tt)<dt1;
        I1=Np(isxf);E1=Ep(isxf);C1=thetap(isxf);
        dm1=2*sqrt(r^2-(sqrt(Y_Au(n_Au)^2+(X_Au(n_Au)+disDec)^2).*sin(C1-tt)).^2);
        
        %Total cross section of the insertion
        Miu1=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(n_Au))*(1+rou_Au(n_Au));
        
        L1=-log(rand(length(I1),1))./Miu1; %Sample path length of particles before interaction
        isxf=L1<dm1; %Ignore path length sampled which is longer than maximum path length in the insertion
        
        I1=I1(isxf);L1=L1(isxf);C1=C1(isxf);dm1=dm1(isxf);E1=E1(isxf);Miu1=Miu1(isxf);
        
        %Calculate the coordinate where scattered photons emitted from
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
        
        %Calculate the number of scattered photons in insertions
        I1=I1.*exp(-exp(interp1(epmma,log(upmma),E1))*1.18.*(l_PMMA-sum(l_Au,2))-sum(l_Au.*miu_Au,2)+l_Au(n_Au).*miu_Au(n_Au))...
            .*(exp(interp1(eli,log(scali),E1))+exp(interp1(eau,log(scaau),E1))*rou_Au(n_Au))./Miu1;
        
        Isca=[Isca;I1];
        C=[C;C1];
        DL=[DL;dl1];
        DM=[DM;dm1];
        LL=[LL;L1];
        Es=[Es;E1];
    end
    clear I1 C1 dl1 dm1 L1 E1
    %--------%
    %Calculate scattered photons in the phantom
    Miu1=exp(interp1(epmma,log(upmma),Ep))*1.18; %Total cross section of the phantom
    
    %Sample path length of particles before interaction
    L1=-log(rand(length(Np),1))./Miu1;
    dl1=zeros(n,1);
    dto=2*(R^2-(disDec*sin(thetap)).^2).^0.5;
    dist=disDec*cos(thetap)-dto/2+L1;
    y_xf=dist.*sin(thetap);
    x_xf=dist.*cos(thetap)-disDec;
    
    %The path length in the insertions should be extracted
    l_Au1=[c_csline(-disDec,0,x_xf,y_xf,X_Au(1),Y_Au(1),r) c_csline(-disDec,0,x_xf,y_xf,X_Au(2),Y_Au(2),r) ...
        c_csline(-disDec,0,x_xf,y_xf,X_Au(3),Y_Au(3),r) c_csline(-disDec,0,x_xf,y_xf,X_Au(4),Y_Au(4),r)];
    l_Au2=zeros(length(L1),4);
    while norm(sum(l_Au1-l_Au2,2)) > 1e-5
        L1=L1+sum(l_Au1-l_Au2,2);
        l_Au2=l_Au1;
        dist=disDec*cos(thetap)-dto/2+L1;
        y_xf=dist.*sin(thetap);
        x_xf=dist.*cos(thetap)-disDec;
        l_Au1=[c_csline(-disDec,0,x_xf,y_xf,X_Au(1),Y_Au(1),r) c_csline(-disDec,0,x_xf,y_xf,X_Au(2),Y_Au(2),r) ...
            c_csline(-disDec,0,x_xf,y_xf,X_Au(3),Y_Au(3),r) c_csline(-disDec,0,x_xf,y_xf,X_Au(4),Y_Au(4),r)];
    end
    
    isxf=L1<dto; %Ignore path length sampled which is longer than maximum path length in the phantom
    L1=L1(isxf);I1=Np(isxf);E1=Ep(isxf);C1=thetap(isxf);l_Au=l_Au1(isxf,:);Miu1=Miu1(isxf);dm1=dto(isxf);dl1=dl1(isxf);
    miu1=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(1))*(1+rou_Au(1));
    miu2=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(2))*(1+rou_Au(2));
    miu3=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(3))*(1+rou_Au(3));
    miu4=(exp(interp1(eli,log(uli),E1))+exp(interp1(eau,log(uau),E1))*rou_Au(4))*(1+rou_Au(4));
    I1=I1.*exp(-sum(l_Au.*[miu1 miu2 miu3 miu4],2))...
        .*exp(interp1(epmma,log(scapmma),E1))*1.18./Miu1;
    
    Isca=[Isca;I1];
    C=[C;C1];
    DL=[DL;dl1];
    DM=[DM;dm1];
    LL=[LL;L1];
    Es=[Es;E1];
    
    dist=disDec*cos(C)+DL-DM/2+LL;
    clear DL DM LL
    
    %Detector 1
    y=dist.*sin(C);
    x=dist.*cos(C)-disDec;
    %Equivalent height of beam
    mh=h*(0.5+0.5*(b+y-keel/2)./a)./(0.5+(b+y-keel/2)./a);
    %Emission angle sampling
    phi1=atan((D/2-x)./(b+y-keel/2));phi2=atan((-x-D/2)./(b+y-keel/2));
    phiflu=phi1+rand(length(Isca),1).*(phi2-phi1);
    %Pixel number
    numdet=floor(((a+b+y).*tan(phiflu)+x)*20+129);
    %Scatter angle
    sita=C+pi/2-phiflu;
    Esca=Es./(1+Es/511.*(1-cos(sita)));
    %Scatter cross section calculated by K-N formula
    alpha=Es/511;
    d_KN=(phi1-phi2).*2*pi.*sin(sita).*1./(1+alpha.*(1-cos(sita))).^2.*((1+cos(sita).^2)/2) ...
        .*(1+alpha.^2.*(1-cos(sita)).^2./((1+cos(sita).^2).*(1+alpha.*(1-cos(sita)))));
    KN=2*pi*((alpha.^3+9*alpha.^2+8*alpha+2)./(alpha.^2.*(1+2*alpha).^2)+ ...
        log(1+2*alpha).*(alpha.^2-2-2*alpha)./(2*alpha.^3));
    persca=d_KN./KN.*mh./(2*pi*(b+y+a)./cos(phiflu));
    Isca1=Isca.*persca;
    
    %Attenuation of scattered photons
    lxf_PMMA=phantom_d(0,-b,x,y,0,0,R);
    lsca_Au=zeros(length(x),4);
    miusca_Au=zeros(length(x),4);
    for n_at=1:4
        lsca_Au(:,n_at)=phantom_d(0,-b,x,y,X_Au(n_at),Y_Au(n_at),r);
        miusca_Au(:,n_at)=(exp(interp1(eli,log(uli),Esca))+exp(interp1(eau,log(uau),Esca))*rou_Au(n_at))*(1+rou_Au(n_at));
    end
    Isca1=Isca1.*exp(-exp(interp1(epmma,log(upmma),Esca))*1.18.*(lxf_PMMA-sum(lsca_Au,2))-sum(miusca_Au.*lsca_Au,2));
    
    %Count by different energy bins
    for n_d=1:256
        ig=Isca1(numdet==n_d);
        eg=Esca(numdet==n_d);
        
        nsca1(1,n_d,nk)=sum(ig((65>eg)&(eg>60)));
        nsca1(2,n_d,nk)=sum(ig((70>eg)&(eg>65)));
        nsca1(3,n_d,nk)=sum(ig((75>eg)&(eg>70)));
        
    end
    
    %Detector 2
    y=-dist.*sin(C);
    x=-(dist.*cos(C)-disDec);
    %Equivalent height of beam
    mh=h*(0.5+0.5*(b+y-keel/2)./a)./(0.5+(b+y-keel/2)./a);
    %Emission angle sampling
    phi1=atan((D/2-x)./(b+y-keel/2));phi2=atan((-x-D/2)./(b+y-keel/2));
    phiflu=phi1+rand(length(Isca),1).*(phi2-phi1);
    %Pixel number
    numdet=floor(((a+b+y).*tan(phiflu)+x)*20+129);
    %Scatter angle
    sita=-C+pi/2+phiflu;
    Esca=Es./(1+Es/511.*(1-cos(sita)));
    %%Scatter cross section calculated by K-N formula
    alpha=Es/511;
    d_KN=(phi1-phi2).*2*pi.*sin(sita).*1./(1+alpha.*(1-cos(sita))).^2.*((1+cos(sita).^2)/2) ...
        .*(1+alpha.^2.*(1-cos(sita)).^2./((1+cos(sita).^2).*(1+alpha.*(1-cos(sita)))));
    KN=2*pi*((alpha.^3+9*alpha.^2+8*alpha+2)./(alpha.^2.*(1+2*alpha).^2)+ ...
        log(1+2*alpha).*(alpha.^2-2-2*alpha)./(2*alpha.^3));
    persca=d_KN./KN.*mh./(2*pi*(b+y+a)./cos(phiflu));
    Isca2=Isca.*persca;
    %Attenuation of scattered photons
    lxf_PMMA=phantom_d(0,-b,x,y,0,0,R);
    lsca_Au=zeros(length(x),4);
    miusca_Au=zeros(length(x),4);
    for n_at=1:4
        lsca_Au(:,n_at)=phantom_d(0,-b,x,y,X_Au(n_at),Y_Au(n_at),r);
        miusca_Au(:,n_at)=(exp(interp1(eli,log(uli),Esca))+exp(interp1(eau,log(uau),Esca))*rou_Au(n_at))*(1+rou_Au(n_at));
    end
    Isca2=Isca2.*exp(-exp(interp1(epmma,log(upmma),Esca))*1.18.*(lxf_PMMA-sum(lsca_Au,2))-sum(miusca_Au.*lsca_Au,2));
    %Count by different energy bins
    for n_d=1:256
        ig=Isca2(numdet==n_d);
        eg=Esca(numdet==n_d);
        
        nsca2(1,n_d,nk)=sum(ig((65>eg)&(eg>60)));
        nsca2(2,n_d,nk)=sum(ig((70>eg)&(eg>65)));
        nsca2(3,n_d,nk)=sum(ig((75>eg)&(eg>70)));
    end
end
%% Total dose calculation
area=2*disDec*tan(phi)*p*(0.5*a+b)/(0.5*a)/(0.2^2);%dose area
N=N.*area*current*t;% total number of incident photons
%Absolute dose conversion
sig1=sum(N(find(E==60):find(E==150)))/sum(Np);
nsca1=nsca1*sig1;
nsca2=nsca2*sig1;
nsca1=random('poisson',nsca1);
nsca2=random('poisson',nsca2);

%Save
save projdata\sca_Au_150keV_double nsca1 nsca2
