clear
clc

set(0,'defaultfigurecolor','w')
load('sca.mat')
load('xrf_Gd_90keV_double.mat')
load('MC.mat')
nscal=nsca1(1,:,1);
nsca=nsca1(2,:,1);
nscar=nsca1(3,:,1);
n_det=n_det1(:,1)';
n1=nscal/sum(nscal+nsca+n_det+nscar)*sum(N1+N2+N3)*1e9;
n2=(nsca+n_det)/sum(nscal+nsca+n_det+nscar)*sum(N1+N2+N3)*1e9;
n3=nscar/sum(nscal+nsca+n_det+nscar)*sum(N1+N2+N3)*1e9;

n1=random('Poisson',n1);
n2=random('Poisson',n2);
n3=random('Poisson',n3);

figure
subplot 131
plot(1:256,n1,1:256,N1*1e9)
legend RMC MCNP5
axis([1 256 0 160])
xlabel({'Pixel';'(a)'})
ylabel('Count')
title('35-40 keV')
subplot 132
plot(1:256,n2,1:256,N2*1e9)
legend RMC MCNP5
axis([1 256 0 160])
xlabel({'Pixel';'(b)'})
ylabel('Count')
title('40-45 keV')
subplot 133
plot(1:256,n3,1:256,N3*1e9)
legend RMC MCNP5
axis([1 256 0 160])
xlabel({'Pixel';'(c)'})
ylabel('Count')
title('45-50 keV')
