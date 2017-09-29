%------------------------%
%MLEM reconstruction algorithm
%------------------------%
clear
clc
set(0,'defaultfigurecolor','w')


l=70;
load('projdata\SMsparse_att_Gd_double.mat')
load('projdata\xrf_Gd_90keV_double.mat')
load('projdata\sca_Gd_90keV_double.mat')
load('projdata\sca_eff_Gd.mat')
n_det=[n_det1;n_det2];
nscal=[squeeze(nsca1(1,:,:));squeeze(nsca2(1,:,:))];
nsca=[squeeze(nsca1(2,:,:));squeeze(nsca2(2,:,:))];
nscar=[squeeze(nsca1(3,:,:));squeeze(nsca2(3,:,:))];


% Greate complete projection data
bin_l=nscal;
bin_m=n_det+nsca;
bin_r=nscar;

% Scatter correction
sca=repmat(sca,2,360);
scal=repmat(scal,2,360);
scar=repmat(scar,2,360);
P=bin_m-(bin_r./scar.*sca+bin_l./scal.*sca)/2;
P(P<0)=0;
imtool(P,[])

%MLEM iteration
p=reshape(P,256*2*360,1);
f=ones(l);
f=reshape(f,l^2,1);
for i=1:15
    f=cij'* (p./(cij*f)).*f./(sum(cij)'); 
    i
end
f=reshape(f,l,l);

imtool(f,[])

