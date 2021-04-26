clear
%% 确定相应参数下E
B=[0:0.01:1]*complex(0,1);
kfasneg=1;eta=40/6;
Etest=[-30:0.1:2];
solverout=Etest*0;
for i=1:length(Etest)
    solverout(i)=dimersolver(Etest(i),kfasneg,eta);
    if solverout(i)>0;
        Eini=fsolve(@(E) dimersolver(E,kfasneg,eta),Etest(i-1));
        break
    end
end
E=0*B+Eini;
kperpc=E;
for i=2:length(E)
    E(i)=fsolve(@(x) dimersolver(x,kfasneg,eta,B(i)),E(i-1));
end
%% 画束缚能随B的变化
figure(1)
subplot(2,1,1)
plot(abs(B),E)
subplot(2,1,2)
plot(abs(B),imag(E))
%% rf谱的计算
omega0=[50:-1:5];
omega1=[5:-0.1:-5];
omega2=[-5:-1:-50];
omega=[omega0,omega1,omega2];


for i=1:length(E)
    for j=1:length(omega)
        
        [a,b,c,d,e]=dimerrf(B(i),E(i),omega(j));
        R(i,j)=a;
        R1(i,j)=d;
        R2(i,j)=e;
        figure(3)
        
        subplot(6,3,1+3*(i-1))
        plot(omega(1:j),R(i,1:j),'linewidth',1.5);
        grid on
        subplot(6,3,2+3*(i-1))
        plot(omega(1:j),R1(i,1:j),'--');
        grid on
        subplot(6,3,3+3*(i-1))
        plot(omega(1:j),R2(i,1:j),'--');
        grid on
       
        
        
       
        drawnow
    end
    
end