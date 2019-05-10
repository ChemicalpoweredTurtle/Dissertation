function [Y]=gauss2(X,Xo,H,w);

% Calculates the gaussian distribution around the center, Xo, for the
% range, X, and half width, w, and the height, H.

[Y]=[H*exp(-1*((X-Xo)/w).^2*(4*log(2)))];

