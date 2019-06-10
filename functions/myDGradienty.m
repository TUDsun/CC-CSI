function Dy = myDGradienty(X, intv)
XR = real(X); XI = imag(X);
method = '2ndOrder';
% DxR = DGradient(XR, intv, 2, '2ndOrder');
DyR = DGradient(XR, intv, 1,method);
% DxI = DGradient(XI, intv, 2, '2ndOrder');
DyI = DGradient(XI, intv, 1,method);
% Dx  = DxR + 1j * DxI;
Dy  = DyR + 1j * DyI;
end
