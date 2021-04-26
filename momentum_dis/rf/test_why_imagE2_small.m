clear
options = optimset('TolFun',1e-10);
B=[0:0.05:1]*complex(0,1);
kfasneg=2;eta=40/6;
%% 找初始值
Etest=[-50:0.001:2];
solverout=Etest*0;
for i=1:length(Etest)
    solverout(i)=dimersolver(Etest(i),kfasneg,eta);
    if solverout(i)>0;
        Eini=fsolve(@(E) dimersolver(E,kfasneg,eta),Etest(i-1));
        break
    end
end
E=0*B+Eini;
%% 迭代，得到每个B对应的2体束缚能
kperpc=E;
for i=2:length(E)
    [e,er]=fsolve(@(x) dimersolver(x,kfasneg,eta,B(i)),E(i-1),options);
    E(i)=e;err(i)=er;
end
%% 束缚能随B的变化
figure(1)
subplot(3,1,1)
plot(abs(B),real(E))
subplot(3,1,2)
plot(abs(B),imag(E))
subplot(3,1,3)
plot(abs(B),abs(err))


%%  固定某个B，算rf谱
i=20;
B=B(i);
E2=E(i);
E=E2;

na=20;
alp=1;
[x1,w1]=lgwt(na,0,1);
x1=flip(x1);w1=flip(w1);
[x2,w2]=lgwt(na,0,1);
k_perp=[x1;1./x2];
wk_perp=[w1;w2./x2.^2];
k_z=k_perp;
wk_z=[w1;w2./x2.^2];

K_perp=repmat(k_perp,1,2*na);
K_z=repmat(k_z',2*na,1);

if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end


Kmod_square=K_perp.^2+K_z.^2;

part1=pi*kfasneg/(eta+1);
x1=B./k_perp./alp;bkp1=-x1-sqrt(x1.^2+1);bkm1=-x1+sqrt(x1.^2+1);

R_kp_up=1./sqrt(1+bkm1.^2);
R_km_up=-1./sqrt(1+bkp1.^2);
R_kp_down=-bkm1./sqrt(1+bkm1.^2);
R_km_down=bkp1./sqrt(1+bkp1.^2);
L_kp_up=-R_km_down;
L_km_up=R_kp_down;
part2temp_1=2*K_perp/(eta+1)./Kmod_square;
part2temp_2=2*R_kp_up.*L_kp_up.*K_perp./(E+Eth-(eta+1)*Kmod_square-2*alp*sqrt(K_perp.^2+B^2)+2*alp*B)+...
    2*R_km_up.*L_km_up.*K_perp./(E+Eth-(eta+1)*Kmod_square+2*alp*sqrt(K_perp.^2+B^2)+2*alp*B);
part2temp=part2temp_1+part2temp_2;
part2=wk_perp'*part2temp*wk_z;
res=part1-part2;

a=2*R_kp_up.*L_kp_up.*K_perp./(E+Eth-(eta+1)*Kmod_square-2*alp*sqrt(K_perp.^2+B^2)+2*alp*B);a=a*wk_z;
b=2*R_km_up.*L_km_up.*K_perp./(E+Eth-(eta+1)*Kmod_square+2*alp*sqrt(K_perp.^2+B^2)+2*alp*B);b=b*wk_z;
%% 画各部分随\k_perp的变化

figure
%+的一支的贡献随k_perp的变化
subplot(3,3,1)
plot(k_perp,real(a)) 
subplot(3,3,2)
plot(k_perp,imag(a))
subplot(3,3,3)
plot(k_perp,abs(a))
%-的一支的贡献
subplot(3,3,4)
plot(k_perp,real(b))
subplot(3,3,5)
plot(k_perp,imag(b))
subplot(3,3,6)
plot(k_perp,abs(b))
%总的rf谱
subplot(3,3,7)
plot(k_perp,real(a+b))
subplot(3,3,8)
plot(k_perp,imag(a+b))
wk_perp'*imag(a+b)
%%