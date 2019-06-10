%% MR_CSI_Brent
% function of the multiplicative-regularized constrast source inversion method


function [Mchi, time, chi, vJ, eTot] = MR_CSI_Brent(Emea, Phi, A, chi, vJ, eIncInv, eTot, pars, epsinv)

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
rchiTV      = pars.rchiTV;
vnu         = cell(1, Nfre);
vg          = cell(1, Nfre);
vgChi       = zeros(size(eIncInv{1}, 1) / NP(ptype), 1);
vnuChi      = cell(1, Nfre);

for ii = 1 : Nfre
    %     eTot{ii}(ri, :) 	= 0;
    vnu{ii}     = zeros(size(eIncInv{1}));
    vg{ii}      = zeros(size(eIncInv{1}));
    vnuChi{ii}  = zeros(size(eIncInv{1}, 1) / NP(ptype), 1);
end

etaS    = cellfun(@(x) norm(x, 'fro'), Emea, ...
    'UniformOutput', false);

etaS    = cellfun(@(x) 1 / (x ^ 2), etaS, ...
    'UniformOutput', false);


vrhotmp = cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
    'UniformOutput', false);

vrho    = cellfun(@(x) x - x, Emea,...
    'UniformOutput', false);



for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end

Rn      = 0;
jj      = 1;
fsb     = 1e3;
tot     = 5;
options = optimset('TolX', 1e-9);
FC      = 1;
nchi    = nnz(eIncInv{1}(:,1));
intrao  = 1;
intv    = 1 / (2 ^ intrao);
logH    = '%5s%10s%10s%10s%13s%13s%13s%13s%13s%13s\n';
logB    = '%5i%10.2e%10.2e%10.2e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e\n';
disp('MRCSI-MF iteration starts \n')
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
    vgs     = cellfun(@(y, z, a) -a * y' * z, Phi, vrho, etaS, ...
        'UniformOutput', false);
    % vgd
    chieTot = cellfun(@(x, y) x * y, Dchi, eTot, ...
        'UniformOutput', false);
    
    chieInc = cellfun(@(x, y) x * y, Dchi, eIncInv, ...
        'UniformOutput', false);
    
    ctmp    = cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
        'UniformOutput', false);
    
    etaD    = cellfun(@(x) x / (x * x + eps), ctmp, ...
        'UniformOutput', false);
    
    vr      = cellfun(@(x, y) x - y, chieTot, vJ, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vr{ii}(ri, :) = 0; end
    
    
    
    
    
    
    
    
    
    
    
    chir    = cellfun(@(x, y) x' * y, Dchi, vr, ...
        'UniformOutput', false);
    
    chiA    = cellfun(@(x, y) x' \ y, A, chir, ...
        'UniformOutput', false); % A'\chir
    
    vgOld   = vg;
    
    vg      = cellfun(@(x, y, z, a) x - a * y + a * z, vgs, vr, chiA, etaD, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vg{ii}(ri, :) = 0; end
    % vNu
    ctmp    = cellfun(@(x) norm(x, 'fro') ^ 2, vgOld, ...
        'UniformOutput', false);
    
    yeta    = cellfun(@(x) x / (x + eps), ctmp, ...
        'UniformOutput', false);
    
    t_v     = cellfun(@(x, y, z, a) x * sum(sum((y - z) .* conj(y))) / (a + eps), yeta, vg, vgOld, ctmp, ...
        'UniformOutput', false); % yeta*sum(sum((vg-vgOld).*conj(vg)))/(ctmp+eps)
    
    vnuOld	= vnu;
    
    vnu     = cellfun(@(x, y, a) x + y * real(a), vg, vnuOld, t_v, ...
        'UniformOutput', false);
    
    eNu     = cellfun(@(x, y) x \ y, A, vnu, ...
        'UniformOutput', false);
    
    %     for ii = 1 : Nfre; eNu{ii}(ri, :) = 0; end
    % va
    chie    = cellfun(@(x, y) x * y, Dchi, eNu, ...
        'UniformOutput', false);
    
    nuchie  = cellfun(@(x, y) x - y, vnu, chie, ...
        'UniformOutput', false);
    
    va      = cellfun(@(x, y, z, a, b, c) -real(sum(conj(x) .* y)) ./ ...
        (z * sum(abs(a * x) .^ 2) + b * sum(abs(c) .^ 2)), ...
        vnu, vg, etaS, Phi, etaD, nuchie, ...
        'UniformOutput', false);
    
    vJ      = cellfun(@(x, y, z) x + y * diag(z), vJ, vnu, va, ...
        'UniformOutput', false);
    
    eTot    = cellfun(@(x, y, z) x + y * diag(z), eTot, eNu, va, ...
        'UniformOutput', false);
    
    
    
    
    
    
    
    chieTot = cellfun(@(x, y) x * y, Dchi, eTot, ...
        'UniformOutput', false);
    
    vrhotmp = cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end
    
    vr      = cellfun(@(x, y) x - y, chieTot, vJ, ...
        'UniformOutput', false);
    
    for ii= 1 : Nfre; vr{ii}(ri, :) = 0; end
    
    
    
    
    
    
    FSC     = cellfun(@(x, a) a * norm(x, 'fro') ^ 2, vrho, etaS, ...
        'UniformOutput', false);
    
    FdC     = cellfun(@(x) norm(x, 'fro') ^ 2, vr, ...
        'UniformOutput', false);
    
    FDC     = cellfun(@(x, y) x * y, etaD, FdC, ...
        'UniformOutput', false);
    
    
    
    
    FS      = sum(cell2mat(FSC)) / Nfre;
    FD      = sum(cell2mat(FDC)) / Nfre;
    %     FC              = sum(cell2mat(FCC)) / Nfre;
    err     = FS + FD;
    
    %% +++++++++++++++++ vgTV +++++++++++++++++
    mchiq       = interp2(mchi{1}, intrao);
    [gx, gy]    = myDGradient(mchiq, intv);
    TVb2        = nchi * (gx .* conj(gx) + gy .* conj(gy) + FD);
    b2chix      = gx ./ (TVb2 + eps);
    b2chiy      = gy ./ (TVb2 + eps);
    bx          = myDGradientx(b2chix, intv);
    by          = myDGradienty(b2chiy, intv);
    vgTVq       = -2 * (bx + by);
    vgTV        = interp2(vgTVq, -intrao);
    vgTV(rchiTV)    = 0;
    
    % +++++++++++++++++ sum(vnuChi, 2) +++++++++++++++++
    vgchiOld    = vgChi;
    
    eTotF       = cellfun(@(x) forward_clb(x, N, Nsrc, ptype), eTot, ...
        'UniformOutput', false);
    
    vrF         = cellfun(@(x) forward_clb(x, N, Nsrc, ptype), vr, ...
        'UniformOutput', false);
    
    
    
    
    
    
    
    
    
    vgChiF      = cellfun(@(x, y, a) 2 * sum(conj(x) .* (a * y), 3).', eTotF, vrF, etaD, ...
        'UniformOutput', false);
    
    vgChiFNor   = cellfun(@(x, y) real(x) + 1j * omega{1} / y * imag(x), vgChiF, omega, ...
        'UniformOutput', false);
    
    vgChiCom    = sum(cell2mat(vgChiFNor), 2) + err * Nfre * vgTV(:);
    
    ETot2       = cellfun(@(x) sum(x .* conj(x), 3).', eTotF, ...
        'UniformOutput', false);
    
    ETot2Nor    = cellfun(@(x, y) x * (omega{1} / y) ^ 2, ETot2, omega, ...
        'UniformOutput', false);
    
    ETot2R      = sum(cell2mat(ETot2), 2);
    ETot2I      = sum(cell2mat(ETot2Nor), 2);
    vgChi       = real(vgChiCom) ./ (ETot2R + eps) + 1j * imag(vgChiCom) ./ (ETot2I + eps); %%%%%%%% normalized !!!!
    vgChi(rchi) = 0;
    stmp        = norm(vgchiOld, 'fro');
    stmp        = stmp ^ 2;
    yeta        = stmp / (stmp + eps);
    t_v         = yeta * sum(sum((vgChi - vgchiOld) .* conj(vgChi))) / (stmp + eps);
    
    vnuChiCom   = vgChi + vnuChi{1} * real(t_v);
    
    vnuChi      = cellfun(@(x) real(vnuChiCom) + 1j * omega{1} / x * imag(vnuChiCom), omega, ...
        'UniformOutput', false);
    
    DnuChi      = cellfun(@(x) back_clb(x, N, ptype), vnuChi, ...
        'UniformOutput', false);
    
    
    % +++++++++++++++++ Step size +++++++++++++++++
    enuchi      = cellfun(@(x, y) x * y, DnuChi, eTot, ...
        'UniformOutput', false);
    
    a00 = sum(cell2mat(FSC));
    a0  = FdC;
    A0  = FDC{1} + FSC{1};
    
    
    
    a1  = cellfun(@(x, y) 2 * real(sum(sum(x .* conj(y)))), enuchi, vr, ...
        'UniformOutput', false);
    
    A1  = cellfun(@(x, y) x * y, etaD, a1, ...
        'UniformOutput', false);
    
    a2  = cellfun(@(x) norm(x, 'fro') ^ 2, enuchi, ...
        'UniformOutput', false);
    
    A2  = cellfun(@(x, y) x * y, etaD, a2, ...
        'UniformOutput', false);
    
    eIncnuchi   = cellfun(@(x, y) x * y, DnuChi, eIncInv, ...
        'UniformOutput', false);
    
    b0  = cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
        'UniformOutput', false);
    
    b1  = cellfun(@(x, y) 2 * real(sum(sum(x .* conj(y)))), eIncnuchi, chieInc, ...
        'UniformOutput', false);
    
    b2  = cellfun(@(x) norm(x, 'fro') ^ 2, eIncnuchi, ...
        'UniformOutput', false);
    
    %     phienuchi       = cellfun(@(x, y) x * y, Phi, enuchi,...
    %     'UniformOutput', false);
    %     c0              = FCC;
    %     c1              = cellfun(@(x, y, z) 2*real(-x * sum(sum(y .* conj(z)))), etaS, phienuchi, vxi,...
    %     'UniformOutput', false);
    %     c2              = cellfun(@(x, y) x * norm(y,'fro')^2, etaS, phienuchi,...
    %     'UniformOutput', false);
    
    d0              = 1;
    vnuChiComq      = interp2(reshape(vnuChiCom, N(1), N(2)), intrao);
    [d_x, d_y]      = myDGradient(vnuChiComq, intv);
    vb1q            = 2 * real((gx .* conj(d_x) + gy .* conj(d_y)) ./ (TVb2 + eps));
    vb1             = interp2(vb1q, -intrao);
    vb1(rchiTV)     = 0;
    d1              = sum(sum(vb1));
    vb2q            = (d_x .* conj(d_x) + d_y .* conj(d_y)) ./ (TVb2 + eps);
    vb2             = interp2(vb2q, -intrao);
    vb2(rchiTV)     = 0;
    d2              = sum(sum(vb2));
    
    c0              = A1{1} * d0 + A0 * d1;
    c1              = 2 * (A1{1} * d1 + A2{1} * d0 + A0 * d2);
    c2              = 3 * (A1{1} * d2 + A2{1} * d1);
    c3              = 4 * A2{1} * d2;
    sol             = roots([c3 c2 c1 c0]);
    sb0             = real(sol(3));
    
    f               = @(x) a00;
    
    for ii = 1 : Nfre
        f = @(x)(a2{ii} * x * x + a1{ii} * x + a0{ii}) / (b2{ii} * x * x + b1{ii} * x + b0{ii}) + f(x);
    end
    
    f               = @(x) f(x) * (d2 * x * x + d1 * x + d0);
    sb              = fminbnd(f, -1e5 * abs(sb0), 1e5 * abs(sb0), options);
    
    % +++++++++++++ Update Chi with Positivity +++++++++++++++++
    chi             = cellfun(@(x, y) x + y * sb, chi, vnuChi, ...
        'UniformOutput', false);
    
    chi             = cellfun(@(x) OrthgPro(x, bgchi), chi, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; chi{ii}(rchi) = 0; end
    
    mchi            = cellfun(@(x) reshape(x, N), chi, ...
        'UniformOutput', false);
    
    Dchi            = cellfun(@(x) back_clb(x, N, ptype), mchi, ...
        'UniformOutput', false);
    
    ChiErr          = norm(epsinv - chi{1}, 'fro') / norm(epsinv, 'fro');
    errOld          = fsb;
    fsb             = f(sb);
    fprintf(logB, jj, FS, FD, FC, fsb, ChiErr, sb0, sb, f(sb), f(sb) - f(sb0))
    
    if mod(jj, interval) == 0
        Mchi{kk}    = mchi{end};
        kk          = kk + 1;
    end
    if (abs(errOld - fsb) <= 1e-5 && terflag)
        tot = tot - 1;
        if tot == 0
            itenum = jj;
        end
    else
        tot = 5;
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





































