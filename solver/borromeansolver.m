function rel=borromeansolver(E,kfasneg,eta,B,alp,cutoff)
% %find mini_eigen_value of the coefficient matirx of borromean binding
%
% define all the things being needed, E is binding energy
% as_lambda_neg is the interaction strength
% alp is the Soc strenght alp=1 means we are in unit of \lambda
% eta is the mass ratio m_a/m_b


% E=real(reE)+complex(0,1)*real(imagE);



switch nargin
    case 3
        B=0;alp=1;invcutoff=0;
    case 4
        alp=1;invcutoff=0;
    case 5
        invcutoff=0;
    case 6
        invcutoff=1/cutoff;
end
    
 


m1=15;
m2=15;
%   define sampling points and weights of q_perp-integral over 0 to 1
[q_perp,wq_perp]=lgwt(m1,0,1);
[q_perpc,wq_perpc]=lgwt(m2,invcutoff,1);
%   define samping points and weights of q_perp-integral over 1 to infty
q_perp_toc=1./q_perpc;
wq_perp_toc=wq_perpc;
%   define the whole sampling points vector and weights vector for q_perp 
Q_perp=[q_perp;q_perp_toc];
Wq_perp=[wq_perp;wq_perp_toc];
%   define sampling points and weights of q_z-integral over 0 to 1
[qz,wqz]=lgwt(m1,0,1);
[qzc,wqzc]=lgwt(m2,0,1);
%   define sampling points and weights of q_z-integral over 1 to infty
qz_toc=1./qzc;
wqz_toc=wqzc;
%   define the whole sampling points vector and weights vector for qz
Qz=[-qz_toc;-qz;qz;qz_toc];
Wqz=[wqz_toc;wqz;wqz;wqz_toc];
%   do the same thing for k-integral
K_perp=Q_perp;
Wk_perp=Wq_perp;
Kz=Qz;
Wkz=Wqz;
intq_perp=[ones(m1,1);q_perp_toc.^2];
intq=[qz_toc.^2;ones(m1,1);ones(m1,1);qz_toc.^2];
[~,~,Intq_perp,Intqz]=ndgrid(ones(m1+m2,1),ones(2*(m1+m2),1),intq_perp,intq);
Int=Intq_perp.*Intqz;
if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end

   
% here f is the first term
f=2*pi*kfasneg./(eta+1);


% because the renormalization term doesn's involve k-integral, seperate
% this term apart from others

% define the mod of k and the mod of q in 4dgrid
[ndK_perp,ndKz,ndQ_perp,ndQz]=ndgrid(K_perp,Kz,Q_perp,Qz);
Kmod=sqrt(ndK_perp.^2+ndKz.^2);
Qmod=sqrt(ndQ_perp.^2+ndQz.^2);

modified_k_perp=sqrt(ndK_perp.^2+B.^2);
modified_q_perp=sqrt(ndQ_perp.^2+B.^2);
x1=B./ndK_perp/alp;bkp1=-x1-sqrt(x1.^2+1);bkm1=-x1+sqrt(x1.^2+1);
x2=B./ndQ_perp/alp;bqp2=-x2-sqrt(x2.^2+1);bqm2=-x2+sqrt(x2.^2+1);
R_kp_up=1./sqrt(1+bkm1.^2);
R_km_up=-1./sqrt(1+bkp1.^2);
R_kp_down=-bkm1./sqrt(1+bkm1.^2);
R_km_down=bkp1./sqrt(1+bkp1.^2);
L_kp_up=-R_km_down;
L_km_up=R_kp_down;

R_qp_up=1./sqrt(1+bqm2.^2);
R_qm_up=-1./sqrt(1+bqp2.^2);
R_qp_down=-bqm2./sqrt(1+bqm2.^2);
R_qm_down=bqp2./sqrt(1+bqp2.^2);
L_qp_up=-R_qm_down;
L_qm_up=R_qp_down;





% the first index is sigma, the second tau: p for plus, m for minus
% define a and b
App=E+2*Eth-(eta+1)*(Kmod.^2+Qmod.^2)-2*alp*(modified_k_perp+modified_q_perp)-2*eta*ndKz.*ndQz+4*alp*B;
Apm=E+2*Eth-(eta+1)*(Kmod.^2+Qmod.^2)-2*alp*(modified_k_perp-modified_q_perp)-2*eta*ndKz.*ndQz+4*alp*B;
Amp=E+2*Eth-(eta+1)*(Kmod.^2+Qmod.^2)-2*alp*(-modified_k_perp+modified_q_perp)-2*eta*ndKz.*ndQz+4*alp*B;
Amm=E+2*Eth-(eta+1)*(Kmod.^2+Qmod.^2)-2*alp*(-modified_k_perp-modified_q_perp)-2*eta*ndKz.*ndQz+4*alp*B;
B=ndK_perp.*ndQ_perp;
B2=B.^2;
X7=2*Int.*ndQ_perp./(eta+1)./(ndQ_perp.^2+ndQz.^2);
% these two are used in the diagnal matrix
X1=Int.*ndQ_perp.*(-2*R_qp_up.*L_qp_up./sqrt(App.^2-4*eta^2*B2) - 2*R_qm_up.*L_qm_up./sqrt((Apm).^2-4*eta^2*B2));
X2=Int.*ndQ_perp.*(-2*R_qp_up.*L_qp_up./sqrt((Amp).^2-4*eta^2*B2) - 2*R_qm_up.*L_qm_up./sqrt((Amm).^2-4*eta^2*B2));
% these four are used to construct the non-diagnal matrix
X3=2*R_kp_up.*L_kp_up.*Int.*ndQ_perp.*(1- App./sqrt((App).^2-4*eta^2*B2))./(2*eta*B);
X4=2*R_kp_up.*L_kp_up.*Int.*ndQ_perp.*(1- Apm./sqrt((Apm).^2-4*eta^2*B2))./(2*eta*B);
X5=2*R_km_up.*L_km_up.*Int.*ndQ_perp.*(1- Amp./sqrt((Amp).^2-4*eta^2*B2))./(2*eta*B);
X6=2*R_km_up.*L_km_up.*Int.*ndQ_perp.*(1- Amm./sqrt((Amm).^2-4*eta^2*B2))./(2*eta*B);
% integrate over qz and q_perp, we have the cutoff term as a number
% integrate X{1} over qz and q_perp, whats left is a 2m*2m matrix
% and put it in a diagnal form
% we have the first coefficent matrix's plus part
Dp=reshape(X1+X7,[],2*(m1+m2))*Wqz;
Dp=reshape(Dp,[],m1+m2)*Wq_perp;
Dp=reshape(Dp,[],1);
Dp=f-Dp;
% here is the minus part
Dm=reshape(X2+X7,[],2*(m1+m2))*Wqz;
Dm=reshape(Dm,[],m1+m2)*Wq_perp;
Dm=reshape(Dm,[],1);
Dm=f-Dm;

D=diag([Dp;Dm]);

[~,~,ndWq_perp,ndWqz]=ndgrid(Wk_perp,Wkz,Wq_perp,Wqz);
Opp=X3.*ndWq_perp.*ndWqz;Opp=reshape(Opp,2*(m1+m2)^2,[]);
Opm=X4.*ndWq_perp.*ndWqz;Opm=reshape(Opm,2*(m1+m2)^2,[]);
Omp=X5.*ndWq_perp.*ndWqz;Omp=reshape(Omp,2*(m1+m2)^2,[]);
Omm=X6.*ndWq_perp.*ndWqz;Omm=reshape(Omm,2*(m1+m2)^2,[]);
O=[Opp Opm;Omp Omm];



M=D+O;

rel=eigs(M,1,'sm');  
end





% something wrong with it: unsolved
% na=10;
% nb=15;
% invcutoff=1/cutoff;
% [x1,w1]=lgwt(na,0,1);
% [x2,w2]=lgwt(nb,0,1);
% kp=[x1;1./x2];
% wkp=[w1;w2];
% kz=[-1./x2;-x1;x1;1./x2];
% wkz=[w2;w1;w1;w2];
% [kp,kz,qp,qz]=ndgrid(kp,kz,kp,kz);
% [~,~,Wqp,Wqz]=ndgrid(wkp,wkz,wkp,wkz);
% intqp=[x1;1./x2.^3];
% intqz=[1./x2.^2;ones(na,1);ones(na,1);1./x2.^2];
% [~,~,Intqp,Intqz]=ndgrid(ones(na+nb,1),ones(2*(na+nb),1),intqp,intqz);
% intMM=Wqp.*Intqp.*Wqz.*Intqz;
% 
% tempM1=kp.^2+kz.^2;
% tempM2=qp.^2+qz.^2;
% denormt=E-(1+eta)*(tempM1+tempM2)-2*eta*kz.*qz;
% App=denormt-kp-qp;
% Apm=denormt-kp+qp;
% Amp=denormt+kp-qp;
% Amm=denormt+kp+qp;
% BM=kp.*qp;
% BM2=BM.^2;
% 
% Dp=2*pi*kfasneg/(1+eta)-(-1./sqrt(App.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM2)).*intMM...
%     -(-1./sqrt(Apm.^2-4*eta^2*BM2)+(1/(1+eta))./tempM2).*intMM;
% Dm=2*pi*kfasneg/(1+eta)-(-1./sqrt(Amp.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM2)).*intMM...
%     -(-1./sqrt(Amm.^2-4*eta^2*BM2)+(1/(1+eta))./tempM2).*intMM;
% Dp=sum(Dp,[3,4]);Dp=reshape(Dp,[],1);
% Dm=sum(Dm,[3,4]);Dm=reshape(Dm,[],1);
% D=diag([Dp;Dm]);
% 
% Opp=((-1-App./sqrt(App.^2-4*eta^2*BM2))./(2*eta*BM)).*intMM;Opp=reshape(Opp,2*(na+nb)^2,[]);
% Opm=((-1-Apm./sqrt(Apm.^2-4*eta^2*BM2))./(2*eta*BM)).*intMM;Opm=reshape(Opm,2*(na+nb)^2,[]);
% Omp=((-1-Amp./sqrt(Amp.^2-4*eta^2*BM2))./(2*eta*BM)).*intMM;Omp=reshape(Omp,2*(na+nb)^2,[]);
% Omm=((-1-Amm./sqrt(Amm.^2-4*eta^2*BM2))./(2*eta*BM)).*intMM;Omm=reshape(Omm,2*(na+nb)^2,[]);
% O=[Opp,Opm;Omp,Omm];
% 
% coeffM=D+O;
% rel=eigs(coeffM,1,'sm');
% 
% end
