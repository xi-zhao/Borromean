clear
B_ini=0;
eta=40/6;
kfasneg_ini=0.5695;
E2_ini=-1.7;
E3_ini=-1.7;
tic
a=find_2b3b_phase_edge(E2_ini,E3_ini,kfasneg_ini,B_ini,eta);
[E2_ini,E3_ini,kfasneg_ini,err_ini]=deal(a(1),a(2),a(3),a(4));
toc
kfasneg=kfasneg_ini-[0:0.001:0.56];
[B,E2,E3,err]=deal(0*kfasneg);
B(1)=B_ini;E2(1)=E2_ini;E3(1)=E3_ini;
tic

for i=2:length(kfasneg)
    a=find_2b3b_phase_edge_scan_kfasneg(E2(i-1),E3(i-1),B(i-1),kfasneg(i),eta);
    [E2(i),E3(i),B(i),err(i)]=deal(a(1),a(2),a(3),a(4));
    figure(1)
    plot(abs(B(1:i)),kfasneg(1:i))
    hold on
    drawnow
    figure(4)
    plot(abs(B(1:i)),err(1:i))
    hold on
    drawnow
    figure(2)
    plot(abs(B(1:i)),real(E2(1:i)))
    hold on
    plot(abs(B(1:i)),real(E3(1:i)))
    drawnow
    figure(3)
    plot(abs(B(1:i)),imag(E2(1:i)))
    hold on
    plot(abs(B(1:i)),imag(E3(1:i)))
    drawnow
    toc
    
end

