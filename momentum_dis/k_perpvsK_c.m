clear

B=[0:0.051:2]*complex(0,1);
k_perp=(1-B.^2).^(1/2);
figure(1)
plot(abs(B),k_perp)
hold on
kfasneg=0.6;eta=40/6;
Etest=[-10:0.1:0];
solverout=Etest*0;
for i=1:length(Etest)
    solverout(i)=borromeansolver(Etest(i),kfasneg,eta);
    if solverout(i)>0;
        Eini=fsolve(@(E) borromeansolver(E,kfasneg,eta),Etest(i-1));
        break
    end
end
E=0*B+Eini;
kperpc=E;
for i=2:length(E)
    E(i)=fsolve(@(x) borromeansolver(x,kfasneg,eta,B(i)),E(i-1));
end
figure(2)
plot(abs(B),E)
for i=1:length(E)
    kperpc(i)=threebodydis(B(i),E(i),kfasneg);
end
figure(1)
hold on
plot(abs(B),kperpc)


