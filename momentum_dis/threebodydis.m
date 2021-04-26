function res=threebodydis(B,E,kfasneg)
%单体哈密顿量的负的那一支的能量的最小值对应的k_perp=(1-B^2)^(1/2)

%需要的参数

eta=40/6;
%拟合的区域
x=[0.2:0.01:1.3];
%得到波函数的整体数据
[a,b]=borromeansolver(E,kfasneg,eta);
%保留负的那一支
eigenvec=b(length(b)/2+1:length(b));
%构造必要的几个系数R，L
invcutoff=0;
[m1,m2]=deal(10,10);
[k_perp,wq_perp]=lgwt(m1,0,1);
[k_perpc,wq_perpc]=lgwt(m2,invcutoff,1);
k_perp=[k_perp;1./k_perpc];

alp=1;
if and(abs(B)<alp,imag(B)==0)
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=0;
end
Kmod_square=k_perp.^2; %这是一列
x1=B./k_perp./alp;bkp1=-x1-sqrt(x1.^2+1);bkm1=-x1+sqrt(x1.^2+1);
R_kp_up=-1./sqrt(1+bkm1.^2);
R_km_up=+1./sqrt(1+bkp1.^2);
R_kp_down=-bkm1./sqrt(1+bkm1.^2);
R_km_down=bkp1./sqrt(1+bkp1.^2);
L_kp_up=-R_km_down;
L_km_up=R_kp_down;

%对波函数重新排序
F_k=reshape(eigenvec,[(m1+m2),2*(m1+m2)]);
F_mk=fliplr(F_k);
%得到未经拟合的波函数
funold=(L_km_up.*L_km_up).^(-1)./(E+2*Eth-2*(Kmod_square-2*...
    sqrt(alp^2*Kmod_square+B^2))).*(R_km_up.*L_km_up).*(F_k+F_mk);
%得到概率分布，并对k_z求和
funold=sum(abs(funold).^2,2)./sum(abs(funold).^2,'all');


% 对k_perp按大小排序
k_perp=[flip(k_perp(1:m1));k_perp(1+m1:m1+m2)];
funold=[flip(funold(1:m1));funold(1+m1:m1+m2)];
% figure
% hold on
% plot(k_perp,funold)
% 开始拟合
dimx=length(x);
xi=k_perp;
%勒让德插值
ximxj=repmat(xi',[length(xi),1])-repmat(xi,[1,length(xi)]);
ximxj=ximxj+eye(length(xi));
ximxjinv=prod(1./ximxj);
xmxj=prod(repmat(x',[1,length(xi)])-repmat(xi',[length(x),1]),2);
fun=((xmxj*ximxjinv)./(repmat(x',[1,length(xi)])-repmat(xi',[dimx,1])))*funold;
fun=sum(fun,2);
[~,ind]=max(fun);
res=x(ind);


figure(3)
hold on

plot(x,fun)
% figure(2)
% phi=-pi:pi/50:pi;
% [K,Phi]=meshgrid(x,phi);
% [Fun,cosphi]=meshgrid(fun,cos(phi));
% Fun=Fun.*cosphi.^2;
% mesh(K.*cos(Phi),K.*sin(Phi),Fun)
end

%

