fn              = length(pars.omega);
va            	= grid3d{fn}.unit.va;
normlz_factor   = pars.omega{fn} * va(PhysQ.omega) * va(PhysQ.eps);

rcT             = pars.rchiTV;
oi              = length(mchi1);

mchi1{oi}(rcT)	= 0;
Epsr           	= real(mchi1{oi});
Sigma          	= -1e3 * imag(mchi1{oi}) * normlz_factor;

figure
subplot(2, 3, 1); myshow2D(yinv, xinv, Epsr, [], fontsize, 'mm')
subplot(2, 3, 4); myshow2D(yinv, xinv, Sigma, [], fontsize, 'mm')


mchi2{oi}(rcT)	= 0;
Epsr           	= real(mchi2{oi});
Sigma          	= -1e3 * imag(mchi2{oi}) * normlz_factor;

% figure
subplot(2, 3, 2); myshow2D(yinv, xinv, Epsr, [], fontsize, 'mm')
subplot(2, 3, 5); myshow2D(yinv, xinv, Sigma, [], fontsize, 'mm')


mchi3{oi}(rcT)	= 0;
Epsr           	= real(mchi3{oi});
Sigma          	= -1e3 * imag(mchi3{oi}) * normlz_factor;

% figure
subplot(2, 3, 3); myshow2D(yinv, xinv, Epsr, [], fontsize, 'mm')
subplot(2, 3, 6); myshow2D(yinv, xinv, Sigma, [], fontsize, 'mm')