clear
B=0;
eta=1;
kfasneg=0.3;
E=linspace(-5,0);
solver=0*E;
for i=1:length(E)
    solver(i)=dimersolver(E(i),kfasneg,eta);
    if solver(i)>0
        break
    end
end
E(i)
figure(1)
plot(E(1:i),solver(1:i))
hold on
plot(E(1:i),solver(1:i)*0)
E0=fsolve(@(x) dimersolver(x,kfasneg,eta),E(i-1));
%%


[k_perp,wk_perp]=lgwt(20,0,4);
[k_z,wk_z]=lgwt(500,-200,200);
[K_perp,K_z]=ndgrid(k_perp,k_z);[Wk_perp,Wk_z]=ndgrid(wk_perp,wk_z);
Ks=K_perp.^2+K_z.^2;
denm=E0-1-eta*Ks-(Ks-2*K_perp);
psi=1./(sqrt(2)*denm);
psivsk_perp=abs(psi).^2*wk_z;

%%
figure
plot(k_perp,psivsk_perp)
figure
mesh(K_perp,K_z,abs(psi).^2)

