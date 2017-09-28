%------------------------%
%本程序为MLEM重建
%------------------------%
clear
clc
set(0,'defaultfigurecolor','w')


l=70;
load('SMsparse_att_Gd_double.mat')
load('xrf_Gd_90keV_double.mat')
load('sca_Gd_90keV_double.mat')
load('sca_eff_90kVp.mat')
n_det=[n_det1;n_det2];
nscal=[squeeze(nsca1(1,:,:));squeeze(nsca2(1,:,:))];
nsca=[squeeze(nsca1(2,:,:));squeeze(nsca2(2,:,:))];
nscar=[squeeze(nsca1(3,:,:));squeeze(nsca2(3,:,:))];




sca=repmat(sca,2,360);
scal=repmat(scal,2,360);
scar=repmat(scar,2,360);

P=n_det+nsca;
Wl=nscal;
Wr=nscar;
P=P-(Wr./scar.*sca+Wl./scal.*sca)/2;
% P=P;
P(P<0)=0;
p=reshape(P,256*2*360,1);
f=ones(l);
f=reshape(f,l^2,1);
for i=1:15
    f=cij'* (p./(cij*f)).*f./(sum(cij)'); 
    i
end
f=reshape(f,l,l);
f=f./3.5054e4*0.4;
% imtool(P,[])
imtool(f,[])

[y x]=meshgrid(1:70);
r1=(x-55.5).^2+(y-35.5).^2<16;
r2=(x-45.5).^2+(y-(35.5+10*sqrt(3))).^2<16;
r3=(x-25.5).^2+(y-(35.5+10*sqrt(3))).^2<16;
r4=(x-15.5).^2+(y-35.5).^2<16;
rb=(x-25.5).^2+(y-(35.5-10*sqrt(3))).^2<16 | (x-45.5).^2+(y-(35.5-10*sqrt(3))).^2<16 ;

recdata=[mean(mean(f(r1))) mean(mean(f(r2))) mean(mean(f(r3))) mean(mean(f(r4))) mean(mean(f(rb)))]
[std(f(r1)) std(f(r2)) std(f(r3)) std(f(r4)) std(f(rb))]
CNR=(recdata-mean(mean(f(rb))))/std(f(rb))

