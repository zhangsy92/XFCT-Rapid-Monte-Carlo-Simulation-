function img=FANFBP_ED_v2(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize)
% img=FANFBP_ED(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize) 
% �Ⱦ������ؽ�
% 2011-3-1 collins
% ��ʽ���Բμ�ҽѧ����ϵͳ���Լ����ִ���
% ��ڲ�����
% p��ͶӰ���ݣ��洢��ʽ��radon�Ľ����ͬ
% theta��ͶӰ�Ƕ�
% DetWidth:̽�������(all)
% DistSourceAxis:Դ����ת���ľ���
% DistSourceDet:Դ��̽��������

% N:�ؽ�ͼ���С
% PixSize:�ؽ����ش�С

% Optimized Huayu Zhang, Jan 2015

ND=size(p,1); %̽��������
NA=size(p,2); %�Ƕ���
if NA~=length(theta) 
    error('����ǶȲ�ƥ��');
end;
DetSize=DetWidth/ND; %̽������Ԫ��С
AngInt=abs(theta(1)-theta(2));%�Ƕȼ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step1��correction
s=((0:ND-1)-0.5*ND+0.5)*DetSize;
cor=DistSourceDet./sqrt(DistSourceDet^2+s.^2);
p=p.*repmat(cor',1,NA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step2��filteration
% 1.use convolution
% h=zeros(2*ND-1,1);
% h(ND)=1/8/DetSize^2;
% h(1:2:2*ND-1)=-1/2/pi^2/DetSize^2./(-ND+1:2:ND-1).^2;
% for i = 1:NA
%    p(:,i) = conv(p(:,i),h,'same');%ע��'same'������
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
p=p(ND:2*ND-1,:);   % ע������������ND������ND+1���;��ͬ�����׷��Ĵ���

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = zeros(N*N,1);
%sin cos tab ע������ϵ
sinbeta = sin(theta+pi/2);
cosbeta = cos(theta+pi/2);

%step3��back projection
%���Է���iradon�Ż������þ���������٣�������Ҫ��չ̽�����Ա��ڴ���̽��������������
y = ((1:N) - 0.5*N - 0.5) * PixSize;
x = ((1:N) - 0.5*N - 0.5) * PixSize;
[X,Y] = meshgrid(x,y);

for t=1:NA
    U=1+(sinbeta(t)*X-cosbeta(t)*Y)/DistSourceAxis;
    c=AngInt./(U.^2);
    s=DistSourceDet*(cosbeta(t)*X+sinbeta(t)*Y)./(DistSourceAxis+sinbeta(t)*X-cosbeta(t)*Y)/DetSize;            
    ind=s+0.5*ND+0.5;%ע��ind�ļ���
    index=floor(ind);
    % remove invalid index
    P_aug = padarray(p, [2,0], 'post'); % for overflow index
    index(index < 1 | ND < index) = ND + 1; 
    lam2=ind-index;lam1=1-lam2; 
    img=img+(lam1(:).*P_aug(index,t)+lam2(:).*P_aug(index+1,t)).*c(:);
end


%������
img=img*DistSourceDet/DistSourceAxis*DetSize;
img = reshape(img,[N,N]);
