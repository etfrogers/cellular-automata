function Ab = function_AbC(dbb,ebb)

 c=3*1e8; lam = 775e-9;   %wavelength 
% Load Ga refractive indices and select at relevant wavelength

% nrload = dlmread('n_poly.txt','\t');krload = dlmread('k_poly.txt','\t');
% nSol = interp1(nrload(:,1),nrload(:,2),lam,'pchip','extrap')+1i*interp1(krload(:,1),krload(:,2),lam,'pchip','extrap');
% nrload2 = dlmread('n_liq.txt','\t');krload2 = dlmread('k_liq.txt','\t');
% nLiq = interp1(nrload2(:,1),nrload2(:,2),lam,'pchip','extrap')+1i*interp1(krload2(:,1),krload2(:,2),lam,'pchip','extrap');
% es = nSol*nSol    
% el =nLiq*nLiq

es=-8.6096 +11.1480i; el=-65.7821 +34.8495i;

degree=0;  f0=c./lam;   w=2*pi*f0;  theta=degree*pi/180;  
da=2.2*1e-9; 
dc=200*1e-9; %da shoule be the same as that used in calcaulation of dynamic process
er0=2.2840; ur0=1;  ern=1;    urn=1;   % parameters of glass substrate  , incident from glass                


layer='ABC';%% A suface melting layer, B composited layer C is the solid layer

para_A=[el 1 da];
para_B=[ebb 1 dbb];
para_C=[es 1 dc];

%==5.0======������ȡ================================================
    para=[para_A;para_B;para_C];
    para_er    = para(:,1);             %erΪ��糣��
    para_ur    = para(:,2);             %urΪ�ŵ���
    para_thick = para(:,3);             %thickΪ���
%==6.0======���Ĥ�еĲ����ֲ�======================================
    for k=1:length(layer)
        type=uint8(layer(k))-64;
        er(k)    = para_er(type);
        ur(k)    = para_ur(type);
        thick(k) = para_thick(type);
    end
    %=======����/������ʵĵ���=====================================
     C_theta0=sqrt(1-sin(theta).^2);
    C_thetan=sqrt(1-sin(theta).^2.*(er0.*ur0)./(ern*urn));
                   %TE��������/������ʵĵ���
        pr0=sqrt(er0)./sqrt(ur0).*C_theta0;
        prn=sqrt(ern)./sqrt(urn).*C_thetan;
   
%==7.0======����������============================================
    M0=eye(2);
    for k=1:length(layer)
        C_theta=sqrt(1-sin(theta).^2.*(er0.*ur0)./(er(k).*ur(k)));
        b=(w/c).*thick(k).*sqrt(er(k)).*sqrt(ur(k)).*C_theta;
                    %TE���Ľ����еĵ���
            pr=sqrt(er(k))./sqrt(ur(k)).*C_theta;
      
        M1=[cos(b) -i./pr.*sin(b); -i.*pr.*sin(b) cos(b)];
        M0=M0*M1;
    end
%==8.0======͸����ϵ������==========================================
    m11=M0(1,1);m12=M0(1,2);m21=M0(2,1);m22=M0(2,2);%��������Ԫ��
    Eout  = 1;
    Ein   = ((m11+m12*prn)+(m21+m22*prn)/pr0)/2*Eout;%���䳡
    Erl   = ((m11+m12*prn)-(m21+m22*prn)/pr0)/2*Eout;%���䳡
    r = Erl/Ein;                    %����ϵ��
    t = sqrt(prn/pr0)*Eout/Ein;     %͸��ϵ��

R  = abs(r)^2;                         %������  
T  = abs(t)^2;    
Ab = 1-R-T ;                        %������

end