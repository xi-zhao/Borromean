clear
warning('off')
options = optimset('TolFun',1e-8,'Display','off');
%找到B=0时的相互作用强度阈值
B_ini=0.546i;
eta=40/6;
kfasneg_ini=0.25748;
E2_ini=-0.070143i;
E3_ini=-1.1214i;

a=find_2b3b_phase_edge(E2_ini,E3_ini,kfasneg_ini,B_ini,eta);
[E2_ini,E3_ini,kfasneg_ini,err_ini]=deal(a(1),a(2),a(3),a(4));


B=B_ini+i*[0:0.00211:0.8];


[kfasneg,E2,E3,err]=deal(0*B);
kfasneg(1)=kfasneg_ini;E2(1)=E2_ini;E3(1)=E3_ini;


for i=2:length(B)
    tic
    a=find_2b3b_phase_edge(E2(i-1),E3(i-1),kfasneg(i-1),B(i),eta);
    [E2(i),E3(i),kfasneg(i),err(i)]=deal(a(1),a(2),a(3),a(4));
    figure(1)
    plot(abs(B(1:i)),kfasneg(1:i))
    hold on
    drawnow
    figure(2)
    plot(abs(B(1:i)),err(1:i))
    hold on
    drawnow
    figure(3)
    plot(abs(B(1:i)),real(E2(1:i)))
    hold on
    plot(abs(B(1:i)),real(E3(1:i)))
    drawnow
    figure(4)
    plot(abs(B(1:i)),imag(E2(1:i)))
    hold on
    plot(abs(B(1:i)),imag(E3(1:i)))
    drawnow
    toc
    
end

