function res=get3bodysol(E3,kfasneg,alp,eta)
% 3body borromean code
% m is the quantum number for z-component of the angular momentum m=1 here
% alp is the rashba soc strength
% kfasneg is the interaction strength
% eta=ma/mb

% method 1 (not possible??)
% ns=10;nb=40;
% [x1,w1]=lgwt(ns,0,1); %k inbetween 0 and 1
% [x2,w2]=lgwt(nb,0,1); %1/k inbetween 1/kc and 1, corresponding to 1<k<kc
% kall=[x1;1./x2];
% kwall=[w1;w2];
% intkp=[x1;1./x2.^3];
% intkz=[ones(ns,1);1./x2.^2];
% 
% kpt=repmat(kall,[1,ns+nb]);
% kzt=repmat(kall',[ns+nb,1]);
% 
% denormt=E3-(1+eta)*(kpt.^2+kz.^2+(kpmat.^2)'+(kzmat.^2)');
% denormpp=denormt-alp*kpmat-alp*kpmat';
% denormpm=denormt-alp*kpmat+alp*kpmat';
% denormmp=denormt+alp*kpmat-alp*kpmat';
% denormmm=denormt+alp*kpmat+alp*kpmat';


%method 2
ns=15;nb=10;
At=1;
[x1,w1]=lgwt(ns,0,At); %k inbetween 0 and 1
[x2,w2]=lgwt(nb,0,1/At); %1/k inbetween 1/kc and 1, corresponding to 1<k<kc
kall=[x1;1./x2]; % k inbetween 0 and kc
kwall=[w1;w2]; % corresponding weight
kzall=[-1./x2;-x1;x1;1./x2]; % kz inbetween -kc and kc
kzwall=[w2;w1;w1;w2]; % kz: corresponding weight
intkp=[x1;1./x2.^3]; % what does this mean? range from 0 to kc^3?
intkz=[1./x2.^2;ones(ns,1);ones(ns,1);1./x2.^2]; % what does this mean? 

kpt=kron(kall,ones(2*(ns+nb),1));  % direct product
kpwt=kron(kwall,ones(2*(ns+nb),1)); % direct product
intkpt=kron(intkp,ones(2*(ns+nb),1)); 
kzt=repmat(kzall,[ns+nb,1]);
kzwt=repmat(kzwall,[ns+nb,1]);
intkzt=repmat(intkz,[ns+nb,1]);
[msize,col]=size(kpt);

kpmat=meshgrid(kpt,ones(1,msize))';
kzmat=meshgrid(kzt,ones(1,msize))';
tempM=meshgrid(kpt.^2+kzt.^2,ones(1,msize))'; % k^2
BM=kpmat.*kpmat';
BM2=BM.^2;

denormt=E3-2*alp^2-(1+eta)*(tempM+tempM')-2*eta*kzmat.*kzmat';
denormpp=denormt-2*alp*kpmat-2*alp*kpmat';
denormpm=denormt-2*alp*kpmat+2*alp*kpmat';
denormmp=denormt+2*alp*kpmat-2*alp*kpmat';
denormmm=denormt+2*alp*kpmat+2*alp*kpmat';

intMM=kpwt.*intkpt.*kzwt.*intkzt;


Ap=2*pi*kfasneg/(1+eta)-(-1./sqrt(denormpp.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM)')*intMM...
     -(-1./sqrt(denormpm.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM)')*intMM;
 Am=2*pi*kfasneg/(1+eta)-(-1./sqrt(denormmp.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM)')*intMM...
     -(-1./sqrt(denormmm.^2-4*eta^2*BM2)+(1/(1+eta))./(tempM)')*intMM;
 Ak=diag([Ap;Am]);
 
 Bpp=((-1-denormpp./sqrt(denormpp.^2-4*eta^2*BM2))./(2*eta*BM))*diag(intMM);
 Bpm=((-1-denormpm./sqrt(denormpm.^2-4*eta^2*BM2))./(2*eta*BM))*diag(intMM);
 Bmp=((-1-denormmp./sqrt(denormmp.^2-4*eta^2*BM2))./(2*eta*BM))*diag(intMM);
 Bmm=((-1-denormmm./sqrt(denormmm.^2-4*eta^2*BM2))./(2*eta*BM))*diag(intMM);
 Bk=[Bpp,Bpm;Bmp,Bmm];
 
 coeffM=Ak+Bk;
 res=eigs(coeffM,1,'sm');

% problematic???
% np=40;nz=40;
% [kp,kpw]=lgwt(np,0,50); %kperp inbetween 0 and kc
% [kz,kzw]=lgwt(nz,-30,30); %kz inbetween 0 and kc (negative values in coefficient)
% 
% % assemble kp and kz vectors for each helicity sector
% kpt=kron(kp,ones(nz,1));
% kpwt=kron(kpw,ones(nz,1));
% kzt=repmat(kz,[np,1]);
% kzwt=repmat(kzw,[np,1]);
% [msize,col]=size(kpt);
% 
% kpmat=meshgrid(kpt,ones(1,msize))';
% kzmat=meshgrid(kzt,ones(1,msize))';
% denormt=E3-(1+eta)*(kpmat.^2+kzmat.^2+(kpmat.^2)'+(kzmat.^2)')-2*eta*kzmat.*kzmat';
% denormpp=denormt-alp*kpmat-alp*kpmat';
% denormpm=denormt-alp*kpmat+alp*kpmat';
% denormmp=denormt+alp*kpmat-alp*kpmat';
% denormmm=denormt+alp*kpmat+alp*kpmat';
% 
% kqmat=kpmat.*kpmat';
% kqmat2=kqmat.^2;
% 
% Ap=2*pi*kfasneg/(1+eta)-(-1./sqrt(denormpp.^2-4*eta^2*kqmat2)+(1/(1+eta))./(kpmat.^2+kzmat.^2)')*(kpt.*kpwt.*kzwt)...
%     -(-1./sqrt(denormpm.^2-4*eta^2*kqmat2)+(1/(1+eta))./(kpmat.^2+kzmat.^2)')*(kpt.*kpwt.*kzwt);
% Am=2*pi*kfasneg/(1+eta)-(-1./sqrt(denormmp.^2-4*eta^2*kqmat2)+(1/(1+eta))./(kpmat.^2+kzmat.^2)')*(kpt.*kpwt.*kzwt)...
%     -(-1./sqrt(denormmm.^2-4*eta^2*kqmat2)+(1/(1+eta))./(kpmat.^2+kzmat.^2)')*(kpt.*kpwt.*kzwt);
% Ak=diag([Ap;Am]);
% 
% Bpp=((-1-denormpp./sqrt(denormpp.^2-4*eta^2*kqmat2))./(2*eta*kqmat))*diag(kpt.*kpwt.*kzwt);
% Bpm=((-1-denormpm./sqrt(denormpm.^2-4*eta^2*kqmat2))./(2*eta*kqmat))*diag(kpt.*kpwt.*kzwt);
% Bmp=((-1-denormmp./sqrt(denormmp.^2-4*eta^2*kqmat2))./(2*eta*kqmat))*diag(kpt.*kpwt.*kzwt);
% Bmm=((-1-denormmm./sqrt(denormmm.^2-4*eta^2*kqmat2))./(2*eta*kqmat))*diag(kpt.*kpwt.*kzwt);
% Bk=[Bpp,Bpm;Bmp,Bmm];
% 
% coeffM=Ak+Bk;
% %coeffM=Ak;
% res=eigs(coeffM,1,'sm');

end

