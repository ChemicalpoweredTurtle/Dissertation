function [Y]=lorzgauss2(X,Xo,H,w,M);

[YG]=gauss2(X,Xo,H,w);
[YL]=lorentz2(X,Xo,H,w);
[Y]=((1-M)*YG)+(M*YL);