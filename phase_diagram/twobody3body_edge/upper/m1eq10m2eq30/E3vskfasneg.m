clear
options = optimset('TolFun',1e-8);
B=0.95*complex(0,1);
tic;
eta=40/6;
kfasneg_ini=0.74911;
er=linspace(-4.3,-4.5,10);
ei=linspace(0,-1,10)*complex(0,1);
[Er,Ei]=ndgrid(er,ei);
E=Er+Ei;
solver=0*Er;
for i=1:length(er)
    for j=1:length(ei)
        solver(i,j)=borromeansolver(E(i,j),kfasneg_ini,eta,B);toc
    end
    
end
figure(1)
mesh(Er,imag(Ei),abs(solver))

kfasneg=linspace(kfasneg_ini,kfasneg_ini-0.8,30);
Ebound=0*kfasneg+E(find(abs(solver)==min(min(abs(solver)))));

error=Ebound;

ab=fsolve(@(x) borromeansolver(x,kfasneg_ini,eta,B),Ebound(1),options);
error(1)=borromeansolver(ab,kfasneg_ini,eta,B);
Ebound(1)=ab;

for j=2:length(kfasneg)
    ab=fsolve(@(x) borromeansolver(x,kfasneg(j),eta,B),Ebound(j-1),options);
    Ebound(j)=ab;
    error(j)=abs(borromeansolver(Ebound(j),kfasneg(j),eta,B));
    figure(2)
    plot(kfasneg(1:j),real(Ebound(1:j)),'b')
    hold on
    drawnow
end
figure(3)
plot(kfasneg(1:j),imag(Ebound(1:j)))
hold on
figure(4)
plot(kfasneg,(error))
hold on
