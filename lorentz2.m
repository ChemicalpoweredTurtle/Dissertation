function [Y]=lorentz2(X,Xo,H,w);

% Calculates the lorentzian distribution around the center, Xo, for the
% range, X, and half width, w, and the height, H.

Y=[H*[4*((X-Xo)/w).^2+1].^-1];