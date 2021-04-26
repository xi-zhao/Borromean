%三体相和两体相的边界
clear
dbstop if error
warning off
eta=40/6;
Bini=0.9;
options=optimset("Display","off");
tic
%找到初始点
Ekfasnegini=fsolve(@(x) E2vsE3(x,eta,Bini),[-0.0668647042433244,0.652052133664724],options);
toc
%扫B的值
[B,E3,kfasneg]=deal(0.9:0.02:1);
B=B*complex(1,0);
E3(1)=Ekfasnegini(1);kfasneg(1)=Ekfasnegini(2);
figure(2)
for i=2:length(B)
    Ekfasneg=fsolve(@(x) E2vsE3(x,eta,B(i)),[E3(i-1),kfasneg(i-1)],options);
    E3(i)=real(Ekfasneg(1));kfasneg(i)=real(Ekfasneg(2));
    toc
    if real(E3(i))>0
        break
    end
    plot((B(1:i)),kfasneg(1:i))
    hold on
    drawnow
end

%画图


%function
function res=E2vsE3(x,eta,B)
res(1)=borromeansolver(real(x(1)),real(x(2)),eta,B);
res(2)=dimersolver(real(x(1)),real(x(2)),eta,B);
end
function res=E2vsE3Imag(x,eta,B)
res(1)= borromeansolver(x(1),real(x(2)),eta,B);
res(2)= dimersolver(x(3),real(x(2)),eta,B);
res(3)=real(x(1))-real(x(3));
end
