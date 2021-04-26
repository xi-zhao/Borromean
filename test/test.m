clear
addpath(genpath('D:\W-WorkBench'))
dbstop if error
kfasneg=0.8;
eta=40/6;
E=linspace(-20,0,10);
out=0*E;
for i=1:length(E)
    out(i)=dimersolvercoeffm(E(i),kfasneg,eta);
end

figure(1)
subplot(3,1,1)
plot(E,real(out))
subplot(3,1,2)
plot(E,imag(out))
subplot(3,1,3)
plot(E,abs(out))