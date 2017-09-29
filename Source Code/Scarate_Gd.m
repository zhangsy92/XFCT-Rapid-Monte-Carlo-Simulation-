%-----------------------------%
% Calculate the theoretical relative number of scattered photons collected 
% by low, middle and high energy bin
%-----------------------------%
clear;clc;
load('data\spec_90kVp_0.3mmCu.mat')
p=0.05;
a=6.35;
n=1e6;
Ep=30+rand(n,1)*60;
Np=interp1(E,N,Ep);
scal=zeros(256,1);
sca=scal;
scar=sca;
for k=1:256
    fprintf('Pixel %d\n',k)
    d=(k-128.5)*p;
    theta=atan(d/a);
    sita=pi/2-theta;
    Esca=Ep./(1+Ep/511.*(1-cos(sita)));
    alpha=Ep/511;
    d_KN=sin(sita).*1./(1+alpha.*(1-cos(sita))).^2.*((1+cos(sita).^2)/2) ...
        .*(1+alpha.^2.*(1-cos(sita)).^2./((1+cos(sita).^2).*(1+alpha.*(1-cos(sita)))));
    Nsca=Np.*d_KN;
    for i=1:n
        if Esca(i)<40 && Esca(i)>35
            scal(k)=scal(k)+Nsca(i);
        end
        if Esca(i)<45 && Esca(i)>40
            sca(k)=sca(k)+Nsca(i);
        end
        if Esca(i)<50 && Esca(i)>45
            scar(k)=scar(k)+Nsca(i);
        end
    end
end

save projdata\sca_eff_Gd scal sca scar