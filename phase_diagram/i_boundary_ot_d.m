%三体相和两体相的边界
clear

dbstop if error
warning off
eta=40/6;

options=optimset("Display","off");
tic
[B,realE,imagE3,imagE2,kfasneg]=deal([[0:0.01:0.41],[0.412:0.002:0.7]]);
B=B*complex(0,1);
%找到初始点(B=0时的E和kfasneg)
a=fsolve(@(x) E2vsE3Imag(x(1),0,0,x(2),eta,B(1)),[-4.2,0.8],options);
realE(1)=a(1);
kfasneg(1)=a(2);
toc
tic
options=optimset('TolFun',1e-6);
error=0*B;
for i=2:10
     
    [a,err]=fsolve(@(x) E2vsE3Imag(real(x(1)),real(x(3)),x(4),real(real(x(2))),eta,B(i)),[realE(i-1),real(imagE3(i-1)),real(imagE2(i-1)),kfasneg(i-1)],options);
    realE(i)=real(a(1));imagE3(i)=real(a(3));imagE2(i)=real(a(4));kfasneg(i)=real(a(2));
    error(i)=err;
    toc
    if  real(realE(i))>0
        break
    end
    figure(5)
    plot(imag(B(1:i)),kfasneg(1:i))
    hold on
    drawnow
    figure(6)
    plot(imag(B(1:i)),realE(1:i))
    hold on 
    drawnow
    figure(7)
    hold on
    plot(imag(B(1:i)),imagE2(1:i))
    drawnow
    figure(8)
    hold on
    plot(imag(B(1:i)),imagE3(1:i))
    drawnow
    figure(9)
    plot(abs(B(1:i)),error(1:i))
end

%画图
% b=0*B;
% c=b;
% for i=1:38
%     res=E2vsE3Imag(realE(i),imagE3(i),imagE2(i),kfasneg(i),eta,B(i));
%     b(i)=res(1);
%     c(i)=res(2);
% end
% figure(1)
% plot(abs(B(1:i)),b(1:i))
% figure(2)
% plot(abs(B(1:i)),c(1:i))


%function

function res=E2vsE3Imag(realE,imagE3,imagE2,kfasneg,eta,B)
a= borromeansolver(real(realE),real(imagE3),real(kfasneg),eta,B);
b= dimersolver(real(realE),real(imagE2),real(kfasneg),eta,B);
res=sqrt(abs(a).^2+abs(b).^2);
end