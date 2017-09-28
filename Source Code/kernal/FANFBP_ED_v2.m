function img=FANFBP_ED_v2(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize)
% img=FANFBP_ED(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize) 
% p£ºProjection data
% theta£ºPrjection angle
% DetWidth:Width of detector (all)
% DistSourceAxis:Distance between source and rotation center
% DistSourceDet:Distance between source and detector

% N:Image size
% PixSize:Pixel size

ND=size(p,1); %Detector size
NA=size(p,2); %projection size
if NA~=length(theta) 
    error('Input angle mismatch');
end;
DetSize=DetWidth/ND; 
AngInt=abs(theta(1)-theta(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step1£ºcorrection
s=((0:ND-1)-0.5*ND+0.5)*DetSize;
cor=DistSourceDet./sqrt(DistSourceDet^2+s.^2);
p=p.*repmat(cor',1,NA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step2£ºfilteration
% 1.use convolution
% h=zeros(2*ND-1,1);
% h(ND)=1/8/DetSize^2;
% h(1:2:2*ND-1)=-1/2/pi^2/DetSize^2./(-ND+1:2:ND-1).^2;
% for i = 1:NA
%    p(:,i) = conv(p(:,i),h,'same');
% end

% 2.use fft
h=zeros(2*ND-1,1);
h(ND)=1/8/DetSize^2;
h(1:2:2*ND-1)=-1/2/pi^2/DetSize^2./(-ND+1:2:ND).^2;
h=[h;zeros(ND-1,1)];
H=fft(h);
p(length(H),1)=0; 
p = fft(p); 
for i = 1:NA
   p(:,i) = p(:,i).*H; 
end
p = ifft(p,'symmetric');     
p=p(ND:2*ND-1,:);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = zeros(N*N,1);
%sin cos tab 
sinbeta = sin(theta+pi/2);
cosbeta = cos(theta+pi/2);

%step3£ºback projection
y = ((1:N) - 0.5*N - 0.5) * PixSize;
x = ((1:N) - 0.5*N - 0.5) * PixSize;
[X,Y] = meshgrid(x,y);

for t=1:NA
    U=1+(sinbeta(t)*X-cosbeta(t)*Y)/DistSourceAxis;
    c=AngInt./(U.^2);
    s=DistSourceDet*(cosbeta(t)*X+sinbeta(t)*Y)./(DistSourceAxis+sinbeta(t)*X-cosbeta(t)*Y)/DetSize;            
    ind=s+0.5*ND+0.5;%
    index=floor(ind);
    % remove invalid index
    P_aug = padarray(p, [2,0], 'post'); % for overflow index
    index(index < 1 | ND < index) = ND + 1; 
    lam2=ind-index;lam1=1-lam2; 
    img=img+(lam1(:).*P_aug(index,t)+lam2(:).*P_aug(index+1,t)).*c(:);
end


%
img=img*DistSourceDet/DistSourceAxis*DetSize;
img = reshape(img,[N,N]);
