fn              = length(pars.omega); 
imageLS       	= reshape(chiinv{fn}, pars.Ninv(1), pars.Ninv(2));
va            	= grid3d{fn}.unit.va;
normlz_factor   = pars.omega{fn} * va(PhysQ.omega) * va(PhysQ.eps);
EpsrLS        	= real(imageLS);
SigmaLS       	= -1e3 * imag(imageLS) * normlz_factor;

dBrange        	= 25;
fontsize       	= 9;
[Xh, Yv]       	= ndgrid(grid3d{1}.l{1 : 2});
yy             	= Xh(:, 1) * grid3d{1}.unitvalue;
xx             	= Yv(1, :).' * grid3d{1}.unitvalue;
yinv           	= 1e3 * yy(pars.ny); 
xinv           	= 1e3 * xx(pars.nx);

clear Xh Yv yy xx

figure
subplot(1, 2, 1); myshow2D(yinv, xinv, EpsrLS, [], fontsize, 'mm')
subplot(1, 2, 2); myshow2D(yinv, xinv, SigmaLS, [], fontsize, 'mm')
plottools('on')