%odinary trimer和borromean相的边界（kfasneg随B的变化）
clear
warning off
dbstop if error
tic
options=optimset("Display","off");
eta=40/6;

%B的取值范围
B=[0:0.0021:1]*complex(0,1);
[kfasneg,eigenvalue]=deal(linspace(0.2,0.4,20));
% 扫kfasneg，确定初始点B（1）处kfasneg的取值
for j = 1:length(kfasneg)
    eigenvalue(j)=dimersolver(0,kfasneg(j),eta,B(1));
    toc
    if eigenvalue(j)>0
        break
    end
end

[E,kfasneg,err]=deal(zeros(1,length(B)));

kfasneg(1)=fsolve(@(x) dimersolver(0,x,eta,B(1)),eigenvalue(j),options);


for i=2:length(E)
    try
        [b,e]=fsolve(@(x) dimersolverErEi(0,real(x(1)),real(x(2)),40/6,B(i)),[E(i-1),kfasneg(i-1)],options);
        E(i)=real(b(1))*complex(0,1);kfasneg(i)=real(b(2));err(i)=e;
        toc
        figure(1)
        plot(abs(B(1:i)),real(E(1:i)))
        drawnow
        figure(2)
        plot(abs(B(1:i)),imag(E(1:i)))
        drawnow
        figure(3)
        plot(abs(B(1:i)),kfasneg(1:i))
        drawnow
        figure(4)
        plot(abs(B(1:i)),err(1:i))
        drawnow
    catch ME
        kfasneg=kfasneg(1:i);
        E=E(1:i)
        B=B(1:i);
        err=err(1:i);
        break;
        
    end
end

