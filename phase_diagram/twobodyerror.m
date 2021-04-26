options=optimset('TolFun',1e-8);

eta=40/6;

for i=1:10
[a,b]=fsolve(@(x) dimersolver(realE(i),imagE2(i),x,eta,B(i)),kfasneg(2),options);
kfasnegerror(i)=kfasneg(i)-a;
error(i)=b;
end
figure(1)
plot(abs(B(1:i)),kfasnegerror)
hold on
figure(2)
plot(abs(B(1:i)),error(1:i));


