clear
options = optimset('TolFun',1e-10);
B=[0:0.1:3]*complex(0,1);
kfasneg=1;eta=40/6;
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
%%
%不同磁场下的能量
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

%%
i=21;
B=B(i);
E2=E(i);
%%
% 测试imag(E2)随kfasneg的变化
kfasnegtest=kfasneg+linspace(0,1);Etest=0*kfasnegtest+E2;
for i=2:length(kfasnegtest)
    Etest(i)=fsolve(@(x) dimersolver(x,kfasnegtest(i),eta,B),Etest(i-1));
end
figure(2)
subplot(1,2,1)
plot(kfasnegtest,real(Etest))
title('Re{E_2} vs (k_f a_s)^{-1}')
hold on
grid on
subplot(1,2,2)
plot(kfasnegtest,imag(Etest))
title('Im{E_2} vs (k_f a_s)^{-1}')
hold on
grid on
%%
alp=1;
cutoff=200;
[k_perp,wk_perp]=lgwt(200,0,20);
[k_z,wk_z]=lgwt(500,-cutoff,cutoff);
[K_perp,K_z]=ndgrid(k_perp,k_z);
for i=1:length(kfasnegtest)
    if and(imag(B)==0,abs(B)<1)
        Eth=real(-alp^2-B^2-2*alp*B);
    elseif  and(abs(B)<alp,imag(B)~=0)
        Eth=real(-alp^2-B^2);
    else
        Eth=real(-4*B);
    end
    Kmod_square=K_perp.^2+K_z.^2;
    
    x1=B./K_perp./alp;bkp1=-x1-sqrt(x1.^2+1);bkm1=-x1+sqrt(x1.^2+1);
    R_kp_up=1./sqrt(1+bkm1.^2);
    R_km_up=-1./sqrt(1+bkp1.^2);
    R_kp_down=-bkm1./sqrt(1+bkm1.^2);
    R_km_down=bkp1./sqrt(1+bkp1.^2);
    L_kp_up=-R_km_down;
    L_km_up=R_kp_down;
    
    funold1=R_km_up./(Etest(i)+Eth-eta*Kmod_square-(Kmod_square-2*B-2*sqrt(alp^2*K_perp.^2+B^2)));
    funold2=R_kp_up./(Etest(i)+Eth-eta*Kmod_square-(Kmod_square-2*B+2*sqrt(alp^2*K_perp.^2+B^2)));
    %+,-基
    fun1=abs(funold1).^2*wk_z;
    fun2=abs(funold2).^2*wk_z;
    %上下基
    funup=funold1.*(-R_km_down)+funold2.*(R_kp_down);
    fundown=funold1.*(R_km_up)+funold2.*(-R_kp_up);
    fun3=abs(funup).^2*wk_z;
    fun4=abs(fundown).^2*wk_z;
    %\up 比 \down
    fracupdown(i)=abs(sum(wk_perp.*k_perp.*fun3)./sum(wk_perp.*k_perp.*fun4));
    fracpm(i)=abs(sum(wk_perp.*k_perp.*fun1)./sum(wk_perp.*k_perp.*fun2));
end
figure(3)
subplot(1,2,1)
plot(kfasnegtest,fracupdown)
hold on
title('\uparrow vs \downarrow')
subplot(1,2,2)
plot(kfasnegtest,fracpm)
hold on
title('+ vs -')


