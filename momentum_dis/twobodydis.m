%% 计算初始B处的E2
clear
options = optimset('TolFun',1e-10);
B=[0:0.001:0.5]*complex(0,1);
kfasneg=0.8;eta=40/6;
Etestr=[-10:0.1:2];
for i=1:length(Etestr)  
    solverout(i)=dimersolver(Etestr(i),kfasneg,eta);  
    if solverout(i)>0
        Eini=(Etestr(i)+Etestr(i-1))/2;
        break
    end
end

figure(1)
plot(Etestr(1:i),solverout(1:i))

E=0*B+Eini;
kperpc=E;
for i=2:length(E)
    [e,er]=fsolve(@(x) dimersolver(x,kfasneg,eta,B(i)),E(i-1),options);
    E(i)=e;err(i)=er;
end
%% 不同磁场下的能量
figure(1)
subplot(3,1,1)
plot(abs(B),real(E))
title('real(E) vs B')
subplot(3,1,2)
plot(abs(B),imag(E))
title('imag(E) vs B')
subplot(3,1,3)
plot(abs(B),abs(err))
title('Err vs B')

%% 设置计算波函数的初始参数 B,E2
i=3;
B=B(i);
E2=E(i);

%% 构造计算波函数的各个部分
alp=1;
if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end
na=30;
invcutoff=0;
[x1,w1]=lgwt(na,0,1);
x1=flip(x1);w1=flip(w1);
[x2,w2]=lgwt(na,invcutoff,1);
k_perp=[x1;1./x2];
wk_perp=[w1;w2./x2.^2];
k_z=k_perp;
wk_z=[w1;w2./x2.^2];
[K_perp,K_z]=ndgrid(k_perp,k_z);
[WK_perp,WK_z]=ndgrid(wk_perp,wk_z);
Kmod_square=K_perp.^2+K_z.^2;

x1=B./k_perp./alp;
bkp1=x1+sqrt(x1.^2+1);bkm1=x1-sqrt(x1.^2+1);
dkp1=x1+conj(sqrt(conj(x1).^2+1));dkm1=x1-conj(sqrt(conj(x1).^2+1));
gamma_kp_up=bkp1./sqrt(1+abs(bkp1).^2);
gamma_km_up=bkm1./sqrt(1+abs(bkm1).^2);
gamma_kp_down=1./sqrt(1+abs(bkp1).^2);
gamma_km_down=1./sqrt(1+abs(bkm1).^2);
bgamma_kp_up=dkp1./sqrt(1+abs(dkp1).^2);
bgamma_km_up=dkm1./sqrt(1+abs(dkm1).^2);
bgamma_kp_down=1./sqrt(1+abs(dkp1).^2);
bgamma_km_down=1./sqrt(1+abs(dkm1).^2);
denomibgamma=bgamma_kp_up.*bgamma_km_down-bgamma_kp_down.*bgamma_km_up;
denomigamma=gamma_kp_up.*gamma_km_down-gamma_kp_down.*gamma_km_up;
L_kp_up=bgamma_km_down./denomibgamma;
L_km_up=-bgamma_kp_down./denomibgamma;
L_kp_down=-bgamma_km_up./denomibgamma;
L_km_down=bgamma_kp_up./denomibgamma;

R_kp_up=gamma_km_down./denomigamma;
R_km_up=-gamma_kp_down./denomigamma;
R_kp_down=-gamma_km_up./denomigamma;
R_km_down=gamma_kp_up./denomigamma;



funold1=R_km_up./(E2+Eth-eta*Kmod_square-(Kmod_square-2*B-2*sqrt(alp^2*K_perp.^2+B^2)));
funold2=R_kp_up./(E2+Eth-eta*Kmod_square-(Kmod_square-2*B+2*sqrt(alp^2*K_perp.^2+B^2)));
%+,-基
C=sum(WK_perp.*K_perp.*(abs(funold1).^2+abs(funold2).^2).*WK_z);
funold1=funold1/sqrt(C);
funold2=funold2/sqrt(C);
fun1=sum(abs(funold1).^2.*WK_z,2);
fun2=sum(abs(funold2).^2.*WK_z,2);
frac=abs(sum(WK_perp.*K_perp.*fun1,'all')./sum(WK_perp.*K_perp.*fun2,'all'))
%上下基
funup=funold1.*(-R_km_down)+funold2.*(R_kp_down);
fundown=funold1.*(R_km_up)+funold2.*(-R_kp_up);
fun3=sum(abs(funup).^2*WK_z,2);
fun4=sum(abs(fundown).^2*WK_z,2);

%\up 比 \down
fracupdown=abs(sum(WK_perp.*K_perp.*fun3,'all')./sum(WK_perp.*K_perp.*fun4,'all'))
fracpm=abs(sum(WK_perp.*K_perp.*fun1,'all')./sum(WK_perp.*K_perp.*fun2,'all'))
%x=0:0.05:2;
% dimx=length(x);
% %generating x-dim
% ximxj=repmat(k_perp',[length(k_perp),1])-repmat(k_perp,[1,length(k_perp)]);
% ximxj=ximxj+eye(length(k_perp));
% ximxjinv=prod(1./ximxj);
% xmxj=prod(repmat(x',[1,length(k_perp)])-repmat(k_perp',[length(x),1]),2);
% fun1=((xmxj*ximxjinv)./(repmat(x',[1,length(k_perp)])-repmat(k_perp',[dimx,1])))*funold1;
% fun1=sum(fun1,2);
% fun2=((xmxj*ximxjinv)./(repmat(x',[1,length(k_perp)])-repmat(k_perp',[dimx,1])))*funold2;
% fun1=sum(fun2,2);

%% 画波函数分量随k_perp的变化
figure(2)
xip=k_perp.^2-2*B+2*sqrt(k_perp.^2+B^2);
xim=k_perp.^2-2*B-2*sqrt(k_perp.^2+B^2);
%波函数
subplot(2,4,[1,5])
plot(k_perp,abs(fun1),'-.b')%-的一支
hold on
plot(k_perp,abs(fun2),'r')%+的一支
title('\psi_{\pm} vs k_{\perp}')
xlim([0 3])
subplot(2,4,[2,6])
plot(k_perp,abs(fun3),'-.r')%\up的一支
hold on
plot(k_perp,abs(fun4),'b')%\down的一支
title('\psi_{\uparrow , \downarrow} vs k_{\perp}')
xlim([0 3])
%\xi_m
subplot(2,4,3)
plot(k_perp,real(xim),'-.b')
title('Re{\xi_{-}} vs k_{\perp}')
xlim([0 3])
subplot(2,4,4)
plot(k_perp,imag(xim),'-.b')
title('Im{\xi_{-}} vs k_{\perp}')
xlim([0 3])
%\xi_p
subplot(2,4,7)
plot(k_perp,real(xip),'r')
title('Re{\xi_{+}} vs k_{\perp}')
xlim([0 3])
subplot(2,4,8)
plot(k_perp,imag(xip),'r')
title('Im{\xi_{+}} vs k_{\perp}')
xlim([0 3])

%% 画各个部分随k_perp的变化
figure(3)
n=7;
subplot(n,3,1)
plot(k_perp,real(R_km_up)*wk_z,'b')
hold on
plot(k_perp,real(R_kp_up)*wk_z,'r')
xlim([0 5])
title('Re{R_{k,\pm}^{\uparrow}} vs k_\perp')

subplot(n,3,2)
plot(k_perp,imag(R_km_up)*wk_z,'b')
hold on
plot(k_perp,imag(R_kp_up)*wk_z,'r')
xlim([0 5])
title('Im{R_{k,\pm}^{\uparrow}} vs k_\perp')

subplot(n,3,3)
plot(k_perp,abs(R_km_up)*wk_z,'b')
hold on
plot(k_perp,abs(R_kp_up)*wk_z,'r')
title('abs{R_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])

subplot(n,3,4)
plot(k_perp,real(L_km_up)*wk_z,'b')
hold on
plot(k_perp,real(L_kp_up)*wk_z,'r')
title('Re{L_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])
subplot(n,3,5)
plot(k_perp,imag(L_km_up)*wk_z,'b')
hold on
plot(k_perp,imag(L_kp_up)*wk_z,'r')
title('Im{L_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])
subplot(n,3,6)
plot(k_perp,abs(L_km_up)*wk_z,'b')
hold on
plot(k_perp,abs(L_kp_up)*wk_z,'r')
title('abs{L_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])
subplot(n,3,7)
plot(k_perp,real(L_km_up.*R_km_up)*wk_z,'b')
hold on
plot(k_perp,real(L_kp_up.*R_kp_up)*wk_z,'r')
title('Re{R_{k,\pm}^{\uparrow}L_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])
subplot(n,3,8)
plot(k_perp,imag(L_km_up.*R_km_up)*wk_z,'b')
hold on
plot(k_perp,imag(L_kp_up.*R_kp_up)*wk_z,'r')
title('Im{R_{k,\pm}^{\uparrow}L_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])
subplot(n,3,9)
plot(k_perp,abs(L_km_up.*R_km_up)*wk_z,'b')
hold on
plot(k_perp,abs(L_km_up.*R_kp_up)*wk_z,'r')
title('abs{R_{k,\pm}^{\uparrow}R_{k,\pm}^{\uparrow}} vs k_\perp')
xlim([0 5])

subplot(n,3,10)
dem=E2+Eth-eta*(K_perp.^2+K_z.^2)-(K_perp.^2+K_z.^2-2*B-2*sqrt(alp^2*K_perp.^2+B^2));
dep=E2+Eth-eta*(K_perp.^2+K_z.^2)-(K_perp.^2+K_z.^2-2*B+2*sqrt(alp^2*K_perp.^2+B^2));
plot(k_perp,real(dem)*wk_z,'b')
hold on
plot(k_perp,real(dep)*wk_z,'r')
title('real part of denominators vs k_\perp')
xlim([0 5])
subplot(n,3,11)
plot(k_perp,imag(dem)*wk_z,'b')
hold on
plot(k_perp,imag(dep)*wk_z,'r')
title('imag part of denominators vs k_\perp')
xlim([0 5])
subplot(n,3,12)
plot(k_perp,abs(dem)*wk_z,'b')
hold on
plot(k_perp,abs(dep)*wk_z,'r')
title('abs of denominators vs k_\perp')
xlim([0 5])
subplot(n,3,13)
dem1=real(E2)+Eth-eta*(K_perp.^2+K_z.^2)-(K_perp.^2+K_z.^2-2*B-2*sqrt(alp^2*K_perp.^2+B^2));
dep1=real(E2)+Eth-eta*(K_perp.^2+K_z.^2)-(K_perp.^2+K_z.^2-2*B+2*sqrt(alp^2*K_perp.^2+B^2));
plot(k_perp,real(dem1)*wk_z,'b')
hold on
plot(k_perp,real(dep1)*wk_z,'r')
title('real part of denominators* vs k_\perp')
xlim([0 5])
subplot(n,3,14)
plot(k_perp,imag(dem1)*wk_z,'b')
hold on
plot(k_perp,imag(dep1)*wk_z,'r')
title('imag part of denominators* vs k_\perp')
xlim([0 5])
subplot(n,3,15)
plot(k_perp,abs(dem1)*wk_z,'b')
hold on
plot(k_perp,abs(dep1)*wk_z,'r')
title('abs of denominators* vs k_\perp')
xlim([0 5])
subplot(n,3,16)
plot(k_perp,real(L_km_up.*R_km_up./dem)*wk_z,'b')
hold on
plot(k_perp,real(L_kp_up.*R_kp_up./dep)*wk_z,'r')
plot(k_perp,real(L_km_up.*R_km_up./dem+L_kp_up.*R_kp_up./dep)*wk_z,'k')
title('real part of the fraction')
xlim([0 5])
subplot(n,3,17)
plot(k_perp,imag(L_km_up.*R_km_up./dem)*wk_z,'b')
hold on
plot(k_perp,imag(L_kp_up.*R_kp_up./dep)*wk_z,'r')
plot(k_perp,imag(L_km_up.*R_km_up./dem+L_kp_up.*R_kp_up./dep)*wk_z,'k')
title('imag part of the fraction')
xlim([0 5])
subplot(n,3,18)
plot(k_perp,abs(L_km_up.*R_km_up./dem)*wk_z,'b')
hold on
plot(k_perp,abs(L_kp_up.*R_kp_up./dep)*wk_z,'r')
plot(k_perp,abs(L_km_up.*R_km_up./dem+L_kp_up.*R_kp_up./dep)*wk_z,'k')
title('abs of the fraction')
xlim([0 5])
subplot(n,3,19)
plot(k_perp,real(L_km_up.*R_km_up./dem1)*wk_z,'b')
hold on
plot(k_perp,real(L_kp_up.*R_kp_up./dep1)*wk_z,'r')
plot(k_perp,real(L_km_up.*R_km_up./dem1+L_kp_up.*R_kp_up./dep1)*wk_z,'k')
title('real part of the fraction*')
xlim([0 5])
subplot(n,3,20)
plot(k_perp,imag(L_km_up.*R_km_up./dem1)*wk_z,'b')
hold on
plot(k_perp,imag(L_kp_up.*R_kp_up./dep1)*wk_z,'r')
plot(k_perp,imag(L_km_up.*R_km_up./dem1+L_kp_up.*R_kp_up./dep1)*wk_z,'k')
title('imag part of the fraction*')
xlim([0 5])
subplot(n,3,21)
plot(k_perp,abs(L_km_up.*R_km_up./dem1)*wk_z,'b')
hold on
plot(k_perp,abs(L_kp_up.*R_kp_up./dep1)*wk_z,'r')
plot(k_perp,abs(L_km_up.*R_km_up./dem1+L_kp_up.*R_kp_up./dep1)*wk_z,'k')
title('abs of the fraction*')
xlim([0 5])
figure(4)
plot(k_perp,imag(L_km_up.*R_km_up./dem+L_kp_up.*R_kp_up./dep)*wk_z,'k')
sum(imag(L_km_up.*R_km_up./dem+L_kp_up.*R_kp_up./dep)*wk_z.*k_perp.*wk_perp)
hold on
plot(k_perp,imag(L_km_up.*R_km_up./dem1+L_kp_up.*R_kp_up./dep1)*wk_z,'.k')
sum(imag(L_km_up.*R_km_up./dem1+L_kp_up.*R_kp_up./dep1)*wk_z.*k_perp.*wk_perp)
title('imag part of the fraction')
xlim([0 5])


%%  计算虚部的加权平均值
xip=K_perp.^2+K_z.^2-2*B+2*sqrt(K_perp.^2+B^2);
xim=K_perp.^2+K_z.^2-2*B-2*sqrt(K_perp.^2+B^2);
imagp=abs(funold1).^2.*imag(xip);
imagm=abs(funold2).^2.*imag(xim);
total_imag=sum(wk_perp.*k_perp.*(imagp+imagm)*wk_z)/C
total_imagp=sum(wk_perp.*k_perp.*(imagp)*wk_z)/C
total_imagm=sum(wk_perp.*k_perp.*(imagm)*wk_z)/C



