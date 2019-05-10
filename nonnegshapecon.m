function [A]=nonnegshapecon(A,Y)
j=1;
    for i=1:Y
        [I J]=max(A(:,i));
        if j==1
            Xo=[J I 10 0.9];
            sigma=[.9 .9 .9 .9];
        else
            Xo=Xo1(i,:);
            sigma=[.9 .9 .9 .9];
        end
        sigmaL=[Xo-sigma.*Xo];
        sigmaU=[Xo+sigma.*Xo];
        options = optimset('Display','off','MaxFunEvals',1E5,'MaxIter',1E5,'TolX',1E-12,'TolFun',1E-12,'LevenbergMarquardt','off','DiffMaxChange',0.1,'DiffMinChange',1E-8,'LargeScale','on'); 
        [Xo1] = lsqnonneg(Xo,[1:length(A)],A(:,i),options);
        [A(:,i)]=peakfit2(Xo1,[1:length(A)]);
        A(:,i)=A(:,i)/norm(A(:,i));
    end
     j=j+1;
end
    
