function img=FANFBP_ED_v2(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize)
% img=FANFBP_ED(p,theta,DetWidth,DistSourceAxis,DistSourceDet,N,PixSize) 
% 等距扇束重建
% 2011-3-1 collins
% 公式可以参见医学成像系统，以及沈乐代码
% 入口参数：
% p：投影数据，存储方式和radon的结果相同
% theta：投影角度
% DetWidth:探测器宽度(all)
% DistSourceAxis:源到旋转中心距离
% DistSourceDet:源到探测器距离

% N:重建图像大小
% PixSize:重建像素大小

% Optimized Huayu Zhang, Jan 2015

ND=size(p,1); %探测器个数
NA=size(p,2); %角度数
if NA~=length(theta) 
    error('输入角度不匹配');
end;
DetSize=DetWidth/ND; %探测器单元大小
AngInt=abs(theta(1)-theta(2));%角度间隔

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step1：correction
s=((0:ND-1)-0.5*ND+0.5)*DetSize;
cor=DistSourceDet./sqrt(DistSourceDet^2+s.^2);
p=p.*repmat(cor',1,NA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step2：filteration
% 1.use convolution
% h=zeros(2*ND-1,1);
% h(ND)=1/8/DetSize^2;
% h(1:2:2*ND-1)=-1/2/pi^2/DetSize^2./(-ND+1:2:ND-1).^2;
% for i = 1:NA
%    p(:,i) = conv(p(:,i),h,'same');%注意'same'的作用
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
p=p(ND:2*ND-1,:);   % 注意这里的起点是ND而不是ND+1，和卷积同样容易犯的错误

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = zeros(N*N,1);
%sin cos tab 注意坐标系
sinbeta = sin(theta+pi/2);
cosbeta = cos(theta+pi/2);

%step3：back projection
%可以仿照iradon优化，利用矩阵运算加速，但是需要扩展探测器以便于处理探测器不够长情形
y = ((1:N) - 0.5*N - 0.5) * PixSize;
x = ((1:N) - 0.5*N - 0.5) * PixSize;
[X,Y] = meshgrid(x,y);

for t=1:NA
    U=1+(sinbeta(t)*X-cosbeta(t)*Y)/DistSourceAxis;
    c=AngInt./(U.^2);
    s=DistSourceDet*(cosbeta(t)*X+sinbeta(t)*Y)./(DistSourceAxis+sinbeta(t)*X-cosbeta(t)*Y)/DetSize;            
    ind=s+0.5*ND+0.5;%注意ind的计算
    index=floor(ind);
    % remove invalid index
    P_aug = padarray(p, [2,0], 'post'); % for overflow index
    index(index < 1 | ND < index) = ND + 1; 
    lam2=ind-index;lam1=1-lam2; 
    img=img+(lam1(:).*P_aug(index,t)+lam2(:).*P_aug(index+1,t)).*c(:);
end


%常数项
img=img*DistSourceDet/DistSourceAxis*DetSize;
img = reshape(img,[N,N]);
