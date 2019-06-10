
function [Mchi, time, chi, vJ, eTot] = CCCSI_BrentMF(...
    Emea, Phi, A, chi, vJ, eIncInv, eTot, pars, epsinv)

interval    = 64;
Mchi        = cell(ceil(pars.itenum / interval) + 1, 1);
kk          = 1;

if isempty(epsinv); epsinv  = chi{1}; end

NP          = 1 : 3;
Nfre        = length(eIncInv);
N           = pars.Ninv;
ri          = pars.rd;
rchi        = pars.rchi;
ptype       = pars.ptype;
% epsinv      = forward_clb(epsinv, N, 1, ptype);
% epsinv      = sum(epsinv, 1) / NP(ptype);
% epsinv      = epsinv.';
bgchi       = pars.bgchi;
omega       = pars.omega;
r           = pars.rmea;
disflag     = pars.disflag;
terflag     = pars.terflag;

mchi        = cellfun(@(x) reshape(x, N), chi, ...
    'UniformOutput', false);

Dchi        = cellfun(@(x) back_clb(x, N, ptype), mchi, ...
    'UniformOutput', false);

Nsrc        = size(Emea{1}, 2);

vnu         = cell(1, Nfre);
vg          = cell(1, Nfre);
vgChi       = zeros(size(eIncInv{1}, 1) / NP(ptype), 1);
vnuChi      = cell(1, Nfre);

for ii = 1 : Nfre
    %     eTot{ii}(ri, :) 	= 0;
    vnu{ii}  	= zeros(size(eIncInv{1}));
    vg{ii}   	= zeros(size(eIncInv{1}));
    vnuChi{ii} 	= zeros(size(eIncInv{1}, 1) / NP(ptype), 1);
end

etaS	= cellfun(@(x) norm(x, 'fro'), Emea, ...
    'UniformOutput', false);

etaS	= cellfun(@(x) 1 / (x ^ 2), etaS, ...
    'UniformOutput', false);


vrhotmp	= cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
    'UniformOutput', false);

vrho 	= cellfun(@(x) x - x, Emea,...
    'UniformOutput', false);

vxi  	= vrho;

for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end

Rn    	= 0;
jj     	= 1;
fsb  	= 1e3;
tot     = 5;
options = optimset('TolX', 1e-9);




logH  	= '%5s%10s%10s%10s%13s%13s%13s%13s%13s%13s\n';
logB  	= '%5i%10.2e%10.2e%10.2e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e\n';
disp('CCCSI-MF iteration starts \n')
fprintf(logH, 'iters', 'FS', 'FD', 'FC', 'err', 'ChiErr', 'sb0', 'sb', 'f', 'f - f0')

%% Iteration starts
if pars.recflag
    diary([pars.str '.txt'])
    diary on
end
startTime = tic;
while jj < pars.itenum + 1
    % +++++++++++++++++ update contrast source +++++++++++++++++
    % vgs
    vgs   	= cellfun(@(y, z, a) -a * y' * z, Phi, vrho, etaS, ...
        'UniformOutput', false);
    % vgd
    chieTot	= cellfun(@(x, y) x * y, Dchi, eTot, ...
        'UniformOutput', false);
    
    chieInc	= cellfun(@(x, y) x * y, Dchi, eIncInv, ...
        'UniformOutput', false);
    
    ctmp  	= cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
        'UniformOutput', false);
    
    etaD   	= cellfun(@(x) x / (x * x + eps), ctmp, ...
        'UniformOutput', false);
    
    vr    	= cellfun(@(x, y) x - y, chieTot, vJ, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vr{ii}(ri, :) = 0; end
    
    vxitmp 	= cellfun(@(x, y, z) x - y * z, Emea, Phi, chieTot, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vxi{ii}(r{ii}) = vxitmp{ii}(r{ii}); end
    
    Phixi 	= cellfun(@(x, y) x' * y, Phi, vxi, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; Phixi{ii}(ri, :) = 0; end
    
    chir 	= cellfun(@(x, y, z, a, b) x' * (a * y - b * z), Dchi, vr, Phixi, etaD, etaS, ...
        'UniformOutput', false); % Dchi'*(etaD*vr-etaS*Phixi)
    
    chiA  	= cellfun(@(x, y) x' \ y, A, chir, ...
        'UniformOutput', false); % A'\chir
    
    vgOld  	= vg;
    
    vg     	= cellfun(@(x, y, z, a) x - a * y + z, vgs, vr, chiA, etaD, ...
        'UniformOutput', false); % vgs - etaD * vr + chiA
    
    for ii = 1 : Nfre; vg{ii}(ri, :) = 0; end
    % vNu
    ctmp   	= cellfun(@(x) norm(x, 'fro') ^ 2, vgOld, ...
        'UniformOutput', false);
    
    yeta  	= cellfun(@(x) x / (x + eps), ctmp, ...
        'UniformOutput', false);
    
    t_v   	= cellfun(@(x, y, z, a) x * sum(sum((y - z) .* conj(y))) / (a + eps), yeta, vg, vgOld, ctmp, ...
        'UniformOutput', false); % yeta*sum(sum((vg-vgOld).*conj(vg)))/(ctmp+eps)
    
    vnuOld 	= vnu;
    
    vnu   	= cellfun(@(x, y, a) x + y * real(a), vg, vnuOld, t_v, ...
        'UniformOutput', false);
    
    eNu  	= cellfun(@(x, y) x \ y, A, vnu, ...
        'UniformOutput', false);
    
    %     for ii = 1 : Nfre; eNu{ii}(ri, :) = 0; end
    % va
    chie  	= cellfun(@(x, y) x * y, Dchi, eNu, ...
        'UniformOutput', false);
    
    nuchie	= cellfun(@(x, y) x - y, vnu, chie, ...
        'UniformOutput', false);
    
    va   	= cellfun(@(x, y, z, a, b, c, d) -real(sum(conj(x) .* y)) ./...
        (c * sum(abs(z * x) .^ 2 + abs(z * a) .^ 2) + d * sum(abs(b) .^ 2)),...
        vnu, vg, Phi, chie, nuchie, etaS, etaD, ...
        'UniformOutput', false);
    
    vJ   	= cellfun(@(x, y, z) x + y * diag(z), vJ, vnu, va, ...
        'UniformOutput', false);
    
    eTot  	= cellfun(@(x, y, z) x + y * diag(z), eTot, eNu, va, ...
        'UniformOutput', false);
    
    
    
    
    
    
    
    chieTot	= cellfun(@(x, y) x * y, Dchi, eTot, ...
        'UniformOutput', false);
    
    vrhotmp	= cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end
    
    vr   	= cellfun(@(x, y) x - y, chieTot, vJ, ...
        'UniformOutput', false);
    
    for ii= 1 : Nfre; vr{ii}(ri, :) = 0; end
    
    vxitmp 	= cellfun(@(x, y, z) x - y * z, Emea, Phi, chieTot, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vxi{ii}(r{ii})  = vxitmp{ii}(r{ii}); end
    
    FSC  	= cellfun(@(x, a) a * norm(x, 'fro') ^ 2, vrho, etaS, ...
        'UniformOutput', false);
    
    FdC  	= cellfun(@(x) norm(x, 'fro') ^ 2, vr, ...
        'UniformOutput', false);
    
    FDC   	= cellfun(@(x, y) x * y, etaD, FdC, ...
        'UniformOutput', false);
    
    FCC   	= cellfun(@(x, y) y * norm(x, 'fro') ^ 2, vxi, etaS, ...
        'UniformOutput', false);
    
    FS   	= sum(cell2mat(FSC)) / Nfre;
    FD    	= sum(cell2mat(FDC)) / Nfre;
    FC    	= sum(cell2mat(FCC)) / Nfre;
    err   	= FS + FD + FC;
    
    %% +++++++++++++++++ vgTV +++++++++++++++++
    %     mchiq           = interp2(mchi{1},intrao);
    %     [gx, gy]        = myDGradient(mchiq,intv);
    %     TVb2            = nchi*(gx.*conj(gx)+gy.*conj(gy)+FD);
    %     b2chix          = gx./(TVb2+eps);
    %     b2chiy          = gy./(TVb2+eps);
    %     bx              = myDGradientx(b2chix,intv);
    %     by              = myDGradienty(b2chiy,intv);
    %     vgTVq           = -2*(bx + by);
    %     vgTV            = interp2(vgTVq,-intrao);
    %     vgTV(rchiTV)    = 0;
    
    % +++++++++++++++++ sum(vnuChi,2) +++++++++++++++++
    vgchiOld 	= vgChi;
    
    eTotF       = cellfun(@(x) forward_clb(x, N, Nsrc, ptype), eTot, ...
        'UniformOutput', false);
    
    vrF         = cellfun(@(x) forward_clb(x, N, Nsrc, ptype), vr, ...
        'UniformOutput', false);
    
    Phixi       = cellfun(@(x, y) x' * y, Phi, vxi, ...
        'UniformOutput', false);
    
    for ii      = 1 : Nfre; Phixi{ii}(ri,:) = 0; end
    
    PhixiF      = cellfun(@(x) forward_clb(x, N, Nsrc, ptype), Phixi, ...
        'UniformOutput', false);
    
    vgChiF      = cellfun(@(x, y, z, a, b) 2 * sum(conj(x) .* (a * y - b * z), 3).', eTotF, vrF, PhixiF, etaD, etaS, ...
        'UniformOutput', false);%%%%%
    
    vgChiFNor	= cellfun(@(x, y) real(x) + 1j * omega{1} / y * imag(x), vgChiF, omega, ...
        'UniformOutput', false);
    
    vgChiCom	= sum(cell2mat(vgChiFNor), 2);
    
    ETot2       = cellfun(@(x) sum(x .* conj(x), 3).', eTotF, ...
        'UniformOutput', false);
    
    ETot2Nor	= cellfun(@(x, y) x * (omega{1} / y) ^ 2, ETot2, omega, ...
        'UniformOutput', false);
    
    ETot2R      = sum(cell2mat(ETot2), 2);
    ETot2I   	= sum(cell2mat(ETot2Nor), 2);
    vgChi    	= real(vgChiCom) ./ (ETot2R + eps) + 1j * imag(vgChiCom) ./ (ETot2I + eps); %%%%%%%% normalized !!!!
    vgChi(rchi)	= 0;
    stmp       	= norm(vgchiOld, 'fro');
    stmp     	= stmp ^ 2;
    yeta      	= stmp / (stmp + eps);
    t_v      	= yeta * sum(sum((vgChi - vgchiOld) .* conj(vgChi))) / (stmp + eps);
    
    vnuChiCom 	= vgChi + vnuChi{1} * real(t_v);
    
    vnuChi   	= cellfun(@(x) real(vnuChiCom) + 1j * omega{1} / x * imag(vnuChiCom), omega, ...
        'UniformOutput', false);
    
    DnuChi  	= cellfun(@(x) back_clb(x, N, ptype), vnuChi, ...
        'UniformOutput', false);
    
    
    % +++++++++++++++++ Step size +++++++++++++++++
    enuchi   	= cellfun(@(x, y) x * y, DnuChi, eTot, ...
        'UniformOutput', false);
    
    a00	= sum(cell2mat(FSC));
    a0	= FdC;
    
    A0	= cellfun(@(x, y) x + y, FDC, FSC, ...
        'UniformOutput', false);
    
    a1	= cellfun(@(x, y) 2 * real(sum(sum(x .* conj(y)))), enuchi, vr, ...
        'UniformOutput', false);
    
    A1	= cellfun(@(x, y) x * y, etaD, a1, ...
        'UniformOutput', false);
    
    a2	= cellfun(@(x) norm(x, 'fro') ^ 2, enuchi, ...
        'UniformOutput', false);
    
    A2	= cellfun(@(x, y) x * y, etaD, a2, ...
        'UniformOutput', false);
    
    eIncnuchi	= cellfun(@(x, y) x * y, DnuChi, eIncInv, ...
        'UniformOutput', false);
    
    b0	= cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
        'UniformOutput', false);
    
    b1	= cellfun(@(x, y) 2 * real(sum(sum(x .* conj(y)))), eIncnuchi, chieInc, ...
        'UniformOutput', false);
    
    b2	= cellfun(@(x) norm(x, 'fro') ^ 2, eIncnuchi, ...
        'UniformOutput', false);
    
    phienuchi	= cellfun(@(x, y) x * y, Phi, enuchi, ...
        'UniformOutput', false);
    
    c0	= cellfun(@(x) norm(x, 'fro') ^ 2, vxi, ...
        'UniformOutput', false);
    
    c1	= cellfun(@(x, y, z) 2 * real(-x * sum(sum(y .* conj(z)))), ...
        etaS, phienuchi,vxi, ...
        'UniformOutput', false);
    
    c2	= cellfun(@(x, y) x * norm(y,'fro') ^ 2, etaS, phienuchi, ...
        'UniformOutput', false);
    
    e0	= sum(cell2mat(A0)) + sum(cell2mat(c0));
    e1	= sum(cell2mat(A1)) + sum(cell2mat(c1));
    e2 	= sum(cell2mat(A2)) + sum(cell2mat(c2));
    sb0	= -e1 / e2 / 2;
    
    
    
    
    
    
    
    
    
    
    f 	= @(x) sum(cell2mat(c2)) * x * x + sum(cell2mat(c1)) * x;
    
    for ii = 1 : Nfre
        f = @(x)(a2{ii} * x * x + a1{ii} * x + a0{ii}) / ...
            (b2{ii} * x * x + b1{ii} * x + b0{ii}) + f(x);
    end
    
    
    sb	= fminbnd(f, -1e5 * abs(sb0), 1e5 * abs(sb0), options);
    
    % +++++++++++++ Update Chi with Positivity +++++++++++++++++
    chi	= cellfun(@(x, y) x + y * sb, chi, vnuChi, ...
        'UniformOutput', false);
    
    chi	= cellfun(@(x) OrthgPro(x, bgchi), chi, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; chi{ii}(rchi) = 0; end
    
    mchi	= cellfun(@(x) reshape(x, N), chi, ...
        'UniformOutput', false);
    
    Dchi	= cellfun(@(x) back_clb(x, N, ptype), mchi, ...
        'UniformOutput', false);
    
    ChiErr	= norm(epsinv - chi{1}, 'fro') / norm(epsinv, 'fro');
    errOld	= fsb;
    fsb  	= f(sb) + a00;
    fprintf(logB, jj, FS, FD, FC, fsb, ChiErr, sb0, sb, f(sb), f(sb) - f(sb0))
    
    if mod(jj, interval) == 0
        Mchi{kk}    = mchi{end};
        kk          = kk + 1;
    end
    if (abs(errOld - fsb) <= 1e-5 && terflag)
        tot	= tot - 1;
        if tot == 0
            pars.itenum = jj;
        end
    else
        tot	= 5;
    end
    if (mod(jj, 32) == 0 && disflag)
        subplot(1, 2, 1); imagesc(real(mchi{end})); axis equal tight; colormap jet; colorbar; drawnow;
        subplot(1, 2, 2); imagesc(-imag(mchi{end})); axis equal tight; colormap jet; colorbar; drawnow;
    end
    jj = jj + 1;
    
end
if pars.recflag; diary off; end
time        = toc(startTime);
Mchi{end}  	= mchi{end};


























