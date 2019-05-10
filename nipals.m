function [P,T]=nipals(X,rank);

% Nonlinear iterative partial least squares
% Renee D. JiJi
% Naval Research Laboratory
% rjiji@ccsalpha3.nrl.navy.mil
% January 17,2001

% Written as described by P. Geladi and B. R. Kowalski
% Analytica Chimica Acta, 185 (1986) 1-17

[I,J]=size(X);
Xnew=X;
t_res=1;
iter=1;
for n=1:rank
   t=Xnew(:,1);
   told=t/norm(t);
   while t_res > 1E-9
      pold=(t'*Xnew)/(t'*t);
      pnew=pold/norm(pold);
      tnew=Xnew*pnew'/(pnew*pnew');
      t_res=abs(sum(told-tnew));
      told=tnew;
   end
   Xold=Xnew;
   t_res=1;
   T(:,n)=tnew;
   P(n,:)=pnew;
   Xnew=Xold-(T(:,n)*P(n,:));   
end

      

