function [sov_fit,subsolv_F,subsolv_Y]=NLLS2(L,H,data,LENGTH,comps);
LENGTH=LENGTH';
[I]=find(LENGTH>L);
range(1)=I(1);
[I]=find(LENGTH>H);
range(2)=I(1);
ydata=data(range(1):range(2),1);
xnew=LENGTH(range(1):range(2),1);
Xo=comps';
sigma=[0.05 .9 .9 1]'*ones(1,(length(comps)/4));
sigma=sigma(:);
sigmaL=[Xo-sigma.*Xo];
sigmaU=[Xo+sigma.*Xo];
N=length(sigma)/4;
sigmaU=reshape(sigmaU,4,N);
sigmaU(4,:)=1;
for i=1:(N)
    if sigmaU(3,i)>50;
    sigmaU(3,i)=50;
    end
end
sigmaU=sigmaU(:);
sigmaL = reshape(sigmaL,4,N);
sigmaL(4,:)=0;
for i=1:(N)
    if sigmaL(3,i)<5;
        sigmaL(3,i)=5;
    end
end
sigmaL=sigmaL(:);
Xo=[Xo;1;-7000];
sigmaL=[sigmaL;-100;-50000];
sigmaU=[sigmaU;100;50000];
options = optimset('MaxFunEvals',1E6,'MaxIter',1E7,'TolX',1E-12,'TolFun',1E-12,'LevenbergMarquardt','off','DiffMaxChange',0.1,'DiffMinChange',1E-8,'LargeScale','on'); 
[sov_fit] = lsqcurvefit(@peakfit3,Xo,xnew,ydata,[sigmaL],[sigmaU],options);
[subsolv_F,subsolv_Y]=peakfit3(sov_fit,LENGTH);
clear Xo sigma sigmaL sigmaU I options range xnew ydata N
