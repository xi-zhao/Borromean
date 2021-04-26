
options = optimset('TolFun',1e-8);
B=0.6*complex(0,1);
eta=40/6;
kfasneg_ini=0.7486;
er=linspace(-3,-6);
ei=linspace(1,-1)*complex(0,1);
[Er,Ei]=ndgrid(er,ei);
E=Er+Ei;
solver=0*Er;
for i=1:length(er)
    for j=1:length(ei)
        solver(i,j)=dimersolver(E(i,j),kfasneg_ini,eta,B);
    end
    
end
figure(1)
mesh(Er,imag(Ei),abs(solver))

kfasneg=linspace(kfasneg_ini,kfasneg_ini-0.5,20);
Ebound=0*kfasneg+E(find(abs(solver)==min(min(abs(solver)))));
error=Ebound;

ab=fsolve(@(x) dimersolver(x,kfasneg_ini,eta,B),Ebound(1),options);
error(1)=dimersolver(ab,kfasneg_ini,eta,B);
Ebound(1)=ab;
figure(2)
for j=2:length(kfasneg)
    ab=fsolve(@(x) dimersolver(x,kfasneg(j),eta,B),Ebound(j-1),options);
    Ebound(j)=ab;
    
    error(j)=abs(dimersolver(Ebound(j),kfasneg(j),eta,B));
    plot(kfasneg(1:j),real(Ebound(1:j)),'k')
    hold on
    drawnow
end
figure(3)
plot(kfasneg(1:j),imag(Ebound(1:j)))
figure(4)
plot(kfasneg,(error))