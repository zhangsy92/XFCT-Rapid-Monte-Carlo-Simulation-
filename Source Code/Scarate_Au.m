%-----------------------------%
% Calculate the theoretical relative number of scattered photons collected 
% by low, middle and high energy bin
%-----------------------------%
clear;clc;
load('data\spec_150kVp_2mmCu.mat')
p=0.05;
a=6.35;
n=1e6;
Ep=60+rand(n,1)*90;
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
        if Esca(i)<65 && Esca(i)>60
            scal(k)=scal(k)+Nsca(i);
        end
        if Esca(i)<70 && Esca(i)>65
            sca(k)=sca(k)+Nsca(i);
        end
        if Esca(i)<75 && Esca(i)>70
            scar(k)=scar(k)+Nsca(i);
        end
    end
end

save projdata\sca_eff_Au scal sca scar