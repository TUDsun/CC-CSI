function Dx = myDGradientx(X,intv)
XR = real(X); XI = imag(X);
method      = '2ndOrder';
DxR = DGradient(XR, intv, 2,method);
% DyR = DGradient(XR, intv, 1, '2ndOrder');
DxI = DGradient(XI, intv, 2,method);
% DyI = DGradient(XI, intv, 1, '2ndOrder');
Dx  = DxR + 1j * DxI;
% Dy  = DyR + 1j * DyI;
end
