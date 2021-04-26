%散射相和borromean相的边界（kfasneg随B的变化）
clear

warning off
dbstop if error
tic
options=optimset("Display","off");
eta=40/6;   
%B的取值范围
B=[0:0.021:1];
B=B*complex(0,1);
kfasneg=zeros(1,length(B));
imagE=kfasneg;imagE(1)=0;
err=B;
kfasneg(1)=fsolve(@(x) borromeansolverErEi(0,0,x,eta,B(1)),0.205,options);
tic
for i=2:length(kfasneg)
    try
        [b,e]=fsolve(@(x) borromeansolverErEi(0,real(x(1)),real(x(2)),40/6,B(i)),[imagE(i-1),kfasneg(i-1)],options);
        kfasneg(i)=real(b(2));imagE(i)=real(b(1));err(i)=e;
        toc
        figure(1)
        plot(abs(B(1:i)),real(kfasneg(1:i)))
        hold on
        drawnow
        figure(2)
        plot(abs(B(1:i)),imagE(1:i))
        hold on
        drawnow
        figure(3)
        plot(abs(B(1:i)),err(1:i))
        drawnow
        toc
    catch ME
        kfasneg=kfasneg(1:i);
        imagE=imagE(1:i)
        B=B(1:i);
        err=err(1:i)
        break;
        
    end
end

