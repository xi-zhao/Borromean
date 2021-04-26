function res=find_2b3b_phase_edge_scan_kfasneg(E2_ini,E3_ini,B_ini,kfasneg,eta)
warning('off')
options = optimset('TolFun',1e-8,'Display','off');

[B,err]=fsolve(@(x) E2vsE3Imag(E2_ini,E3_ini,kfasneg,eta,x),B_ini,options);

E3=fsolve(@(x) borromeansolver(x,kfasneg,eta,B),E3_ini,options);
E2=fsolve(@(x) dimersolver(x,kfasneg,eta,B),E2_ini,options);

res(1)=E2;
res(2)=E3;
res(3)=B;
res(4)=err;

    function res=E2vsE3Imag(E2,E3,kfasneg,eta,B)
        options = optimset('TolFun',1e-8,'Display','off');
        Eb3=fsolve(@(x) borromeansolver(x,kfasneg,eta,B),E3,options);
        Eb2=fsolve(@(x) dimersolver(x,kfasneg,eta,B),E2,options);
        res=real(Eb2-Eb3);
       
    end
end





