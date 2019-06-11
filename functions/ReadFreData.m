function [dat, Phi, Ainv, Einc, grid3d, pars, pars_LSM] ...
    = ReadFreData(fre, rawdat, invdom, regSize, para)

pars_LSM        = [];
rawdat          = rawdat(rawdat(:, 3) == fre, :);
m_unit          = para.m_unit;
r               = 1 : size(rawdat, 1);
NRX             = size(rawdat, 1) / para.NTX;
r               = reshape(r, NRX, para.NTX);
% r(1 : 2 : end, :) = [];
r               = [];
r               = r(:);
rawdat(r, :)    = [];
NRX             = size(rawdat, 1) / para.NTX;
calind          = floor(NRX/2) + 1;
Rindex          = reshape(rawdat(:, 2), NRX, para.NTX);
Runiq           = unique(Rindex);
[~, Rindex]   	= ismember(Rindex, Runiq);

Ftot                = rawdat(:, 4) + 1j * rawdat(:, 5);
Finc                = rawdat(:, 6) + 1j * rawdat(:, 7);
rawdat(:, 3)        = Ftot - Finc;
rawdat(:, 4 : end)	= [];

% Trad            = deg2rad(para.TxInterval * (rawdat(:, 1) - 1) - 2.5);
% Rrad            = deg2rad(5 * (rawdat(:, 2) - 1));
% Rraduniq        = deg2rad(5 * (Runiq - 1));
Trad            = deg2rad(rawdat(:, 1));
Rrad            = deg2rad(rawdat(:, 2));
Rraduniq        = deg2rad(Runiq);

rawdat(:, 1 : 2)	= [];

if para.pt == PT.TE
    rawdat	= [rawdat .* cos(Rrad + pi / 2), rawdat .* sin(Rrad + pi / 2)].';
    Ftot    = [Ftot .* cos(Rrad + pi / 2), Ftot .* sin(Rrad + pi / 2)].';
    FInc    = Finc;
    Finc    = [FInc .* cos(Rrad + pi / 2), FInc .* sin(Rrad + pi / 2)].';
    dat     = reshape(rawdat(:), 2 * NRX, para.NTX);
    Rindex  = [2 * Rindex(:) - 1, 2 * Rindex(:)].';
    Rindex  = reshape(Rindex(:), numel(Rindex) / para.NTX, para.NTX);
    Ftot    = reshape(Ftot(:), 2 * NRX, para.NTX);
    Finc    = reshape(Finc(:), 2 * NRX, para.NTX);
    FInc    = reshape(FInc(:), NRX, para.NTX);
else
    dat     = reshape(rawdat, NRX, para.NTX);
    Ftot    = reshape(Ftot(:), NRX, para.NTX);
    Finc    = reshape(Finc(:), NRX, para.NTX);
end
Tx              = para.Tr * cos(Trad);
Ty              = para.Tr * sin(Trad);
Rx              = para.Rr * cos(Rrad);
Ry              = para.Rr * sin(Rrad);
Rxuniq          = para.Rr * cos(Rraduniq) / m_unit;
Ryuniq          = para.Rr * sin(Rraduniq) / m_unit;
ind             = (1 : NRX) + 0 * NRX;
if para.showconfg
    figure
    plot(Tx(1), Ty(1), 'kd','MarkerFaceColor','k')
    hold on
    plot(Rx(ind), Ry(ind), 'ko')
    axis equal
    hold on
    plot(Tx(NRX+1), Ty(NRX+1), 'kd')
    xlabel('$x$/m', 'interpreter', 'latex')
    ylabel('$y$/m', 'interpreter', 'latex')
    lgnd = legend('Current Transmitter', 'Current Receiver', 'Next Receiver');
    set(lgnd, 'color', 'white', 'Location', 'best', 'FontName', 'Times New Roman');
    grid on
end

%% Create Ms and A
lambda          = PhysC.c0 / (fre * 1e9);
wvlenn          = lambda / m_unit;
regGrid         = round(regSize / m_unit);
d_pml           = 5;
epsr            = 1;
regGrid         = [regGrid(:, 1) - d_pml, regGrid(:, 2) + d_pml];

if para.pt == PT.TE
    MsAname     = ['MsA' num2str(fre) 'GHzTE.mat'];
    Phiname     = ['Phi' num2str(fre) 'GHzTE.mat'];
else
    MsAname     = ['MsA' num2str(fre) 'GHz.mat'];
    Phiname     = ['Phi' num2str(fre) 'GHz.mat'];
end

if exist(Phiname, 'file')
    load(Phiname);
else
    posRxJ          = round(Rxuniq);
    posRyJ          = round(Ryuniq);
    posRxM          = round(Rxuniq - 0.5) + 0.5;
    posRyM          = round(Ryuniq - 0.5) + 0.5;

    switch para.pt
        case PT.FULL
            Src_J           = cell(length(Runiq), 1);
            Src_M           = cell(length(Runiq), 1);
            for ii = 1 : length(Runiq)
                Src_J{ii}   =  PointSrc(Axis.z, [posRxJ(ii), posRyJ(ii), 0.5]);
                Src_M{ii}	=  PointSrc(Axis.z, [posRxM(ii), posRyM(ii),   0]);
            end
            SolverOpts.returnAandb  = false;
            [osc, grid3d, ~, ~, ~] = FDFDpars(...
                'OSC', m_unit, wvlenn, ...
                'DOM', {'vacuum', 'none', epsr}, [regGrid; 0, 1], 1, BC.p, [d_pml d_pml 0], ...
                'SRCJ', Src_J{:}, 'SRCM', Src_M{:}, SolverOpts, para.pt);
        case PT.TE
            Src_M = cell(length(Runiq), 1);
            SolverOpts.returnAandb  = false;
            [osc, grid3d, ~, ~, ~] = FDFDpars(...
                'OSC', m_unit, wvlenn, ...
                'DOM', {'vacuum', 'none', epsr}, [regGrid; 0, 1], 1, BC.p, [d_pml d_pml 0], ...
                'SRCM', Src_M{:}, SolverOpts, para.pt);
        case PT.TM
            Src_J = cell(length(Runiq), 1);
            SolverOpts.returnAandb  = false;
            [osc, grid3d, ~, ~, ~] = FDFDpars(...
                'OSC', m_unit, wvlenn, ...
                'DOM', {'vacuum', 'none', epsr}, [regGrid; 0, 1], 1, BC.p, [d_pml d_pml 0], ...
                'SRCJ', Src_J{:}, SolverOpts, para.pt);
    end

    xlchi           = round(invdom(1) / m_unit);
    xhchi           = round(invdom(2) / m_unit);
    ylchi           = round(invdom(3) / m_unit);
    yhchi           = round(invdom(4) / m_unit);
    zlchi           = round(invdom(5) / m_unit);
    zhchi           = round(invdom(6) / m_unit);
    
    xpuni           = cell2mat(grid3d.l(1, 1));
    xduni           = cell2mat(grid3d.l(1, 2));
    ypuni           = cell2mat(grid3d.l(2, 1));
    yduni           = cell2mat(grid3d.l(2, 2));
    zpuni           = cell2mat(grid3d.l(3, 1));
    zduni           = cell2mat(grid3d.l(3, 2));
    
    [~, ixlp]       = min(abs(xlchi - xpuni));
    [~, ixhp]       = min(abs(xhchi - xpuni));
    [~, iylp]       = min(abs(ylchi - ypuni));
    [~, iyhp]       = min(abs(yhchi - ypuni));
    [~, izlp]       = min(abs(zlchi - zpuni));
    [~, izhp]       = min(abs(zhchi - zpuni));
    
    %% Compute A
    dl              = xpuni(ixlp) - xpuni(ixlp - 1);
    dlvac           = 1;
    dlbg            = 1;
    if(abs(dl - dlvac) < abs(dl - dlbg))
        dl = dlvac;
    else
        dl = dlbg;
    end
    if abs(xpuni(ixlp) - dl - xpuni(ixlp - 1)) < 1e-3 * dl % xpuni(ixlp) - dl == xpuni(ixlp - 1)
        xlinv = xpuni(ixlp);
    else
        xlinv = xpuni(ixlp) - dl;
        ixlp  = ixlp - 1;
    end
    
    dl = xpuni(ixhp + 1) - xpuni(ixhp);
    if(abs(dl - dlvac) < abs(dl - dlbg))
        dl = dlvac;
    else
        dl = dlbg;
    end
    if abs(xpuni(ixhp) + dl - xpuni(ixhp + 1)) < 1e-3 * dl % xpuni(ixhp) + dl == xpuni(ixhp + 1)
        xhinv = xpuni(ixhp);
    else
        xhinv = xpuni(ixhp) + dl;
        ixhp  = ixhp + 1;
    end
    
    dl = ypuni(iylp) - ypuni(iylp - 1);
    if(abs(dl - dlvac) < abs(dl - dlbg))
        dl = dlvac;
    else
        dl = dlbg;
    end
    if abs(ypuni(iylp) - dl - ypuni(iylp - 1)) < 1e-3 * dl % ypuni(iylp) - dl == ypuni(iylp - 1)
        ylinv = ypuni(iylp);
    else
        ylinv = ypuni(iylp) - dl;
        iylp  = iylp - 1;
    end
    
    dl = ypuni(iyhp + 1) - ypuni(iyhp);
    if(abs(dl - dlvac) < abs(dl - dlbg))
        dl = dlvac;
    else
        dl = dlbg;
    end
    if abs(ypuni(iyhp) + dl - ypuni(iyhp + 1)) < 1e-3 * dl % ypuni(iyhp) + dl == ypuni(iyhp + 1)
        yhinv = ypuni(iyhp);
    else
        yhinv = ypuni(iyhp) + dl;
        iyhp  = iyhp + 1;
    end
    
    %     dl = zpuni(izlp)-zpuni(izlp-1);
    %     if(abs(dl-dlvac)<abs(dl-dlbg))
    %         dl = dlvac;
    %     else
    %         dl = dlbg;
    %     end
    %     if zpuni(izlp)-dl == zpuni(izlp-1)
    %         zlinv = zpuni(izlp);
    %     else
    %         zlinv = zpuni(izlp)-dl;
    %         izlp  = izlp-1;
    %     end
    %
    %     dl = zpuni(izhp+1)-zpuni(izhp);
    %     if(abs(dl-dlvac)<abs(dl-dlbg))
    %         dl = dlvac;
    %     else
    %         dl = dlbg;
    %     end
    %     if zpuni(izhp)+dl == zpuni(izhp+1)
    %         zhinv = zpuni(izhp);
    %     else
    %         zhinv = zpuni(izhp)+dl;
    %         izhp  = izhp+1;
    %     end
    %%
    vac     = Material('vacuum', 'non', epsr);
    domain  = Domain([...
        xlinv, xhinv; ...
        ylinv, yhinv; ...
        0, 1], dl);
    % srcx        = round((xlinv + xhinv) / 2);
    % srcy        = round((ylinv + yhinv) / 2);
    srcx    = xpuni(round((ixlp + ixhp) / 2));
    srcy    = ypuni(round((iylp + iyhp) / 2));
    srcz    = 0.5;

    DomVacInv               = EMObject(domain, vac);
    SolverOpts.returnAandb  = true;

    [~, grid3dinv, ~, Ainv, ~] = FDFDpars(...
        'OSC', m_unit, wvlenn, ...
        'DOM', DomVacInv, BC.p, [d_pml d_pml 0], ...
        'SRCJ', PointSrc(Axis.z, [srcx, srcy, srcz]), ...
        SolverOpts, para.pt);
    
    %% Compute Phi
    nx              = ixlp : ixhp - 1;
    ny              = iylp : iyhp - 1;
    nz              = izlp : izhp;
    Ninv            = [length(nx), length(ny), length(nz)];
    prodN           = prod(Ninv);
    prodN2          = prodN * 2;
    pars.nx         = nx;
    pars.ny         = ny;
    pars.nz         = nz;
    pars.Ninv       = Ninv;
    omega           = osc.in_omega0();
    pars.omega      = omega;
    k               = omega * sqrt(1);

















    switch para.pt
        case PT.FULL
            
        case PT.TE
            Phi             = zeros(length(Rxuniq) * 2, prodN2);
            [NY1, NX1, ~] 	= meshgrid(ypuni(ny), xduni(nx), zpuni(nz));
            Dx1             = repmat(Rxuniq, 1, numel(NX1)) - repmat(NX1(:).', numel(Rxuniq), 1);
            Dy1             = repmat(Ryuniq, 1, numel(NY1)) - repmat(NY1(:).', numel(Ryuniq), 1);
            [NY2, NX2, ~] 	= meshgrid(yduni(ny), xpuni(nx), zpuni(nz));
            Dx2             = repmat(Rxuniq, 1, numel(NX2)) - repmat(NX2(:).', numel(Rxuniq), 1);
            Dy2             = repmat(Ryuniq, 1, numel(NY2)) - repmat(NY2(:).', numel(Ryuniq), 1);
            
            Dxs1            = Dx1 .^ 2;
            Dys1            = Dy1 .^ 2;
            Dxs2            = Dx2 .^ 2;
            Dys2            = Dy2 .^ 2;
            
            Rs1             = Dxs1 + Dys1;
            Rs2             = Dxs2 + Dys2;
            
            R1              = sqrt(Rs1);
            R2              = sqrt(Rs2);
            
            C11             = -k / (omega * 4) ./ (R1 + eps);
            C12             = k ^ 2 / (omega * 4) .* Dx2 .* Dy2 ./ (Rs2 + eps);
            C21             = k ^ 2 / (omega * 4) .* Dx1 .* Dy1 ./ (Rs1 + eps);
            C22             = -k / (omega * 4) ./ (R1 + eps);
            
            Phi(1 : 2 : end, 1 : 2 : end) = C11 .* (besselh(1, 1, -k * R1) + k * besselh(2, 1, -k * R1) .* Dys1 ./ (R1 + eps));
            Phi(1 : 2 : end, 2 : 2 : end) = C12 .* besselh(2, 1, -k * R2);
            Phi(2 : 2 : end, 1 : 2 : end) = C21 .* besselh(2, 1, -k * R1);
            Phi(2 : 2 : end, 2 : 2 : end) = C22 .* (besselh(1, 1, -k * R2) + k * besselh(2, 1, -k * R2) .* Dxs2 ./ (R2 + eps));
            
            Phi             = Phi * omega^2 / (-1j * omega);
            
            k               = k / m_unit;
            Einc            = zeros(prodN2, para.NTX);
            
            for ii = 1 : para.NTX
                Txiic           = Tx(NRX * ii - NRX + 1, :);
                Tyiic           = Ty(NRX * ii - NRX + 1, :);
                TyiiC           = repmat(Tyiic, size(NX1));
                TxiiC           = repmat(Txiic, size(NX1));
                yyRC            = NY1 * m_unit - TyiiC;
                xxRC            = NX1 * m_unit - TxiiC;
                RRC             = sqrt(yyRC .^ 2 + xxRC .^ 2);
                
                EPhix           = -1j * k / 4 * yyRC .* besselh(1, 1, -k * RRC) ./ (RRC + eps);
                
                TyiiC           = repmat(Tyiic, size(NX2));
                TxiiC           = repmat(Txiic, size(NX2));
                yyRC            = NY2 * m_unit - TyiiC;
                xxRC            = NX2 * m_unit - TxiiC;
                RRC             = sqrt(yyRC .^ 2 + xxRC .^ 2);
                
                EPhiy           = +1j * k / 4 * xxRC .* besselh(1, 1, -k * RRC) ./ (RRC + eps);
                
                EPhi            = [EPhix(:).'; EPhiy(:).'];
                EPhi            = EPhi(:);
                
                ind             = (1 : NRX) + (ii - 1) * NRX;
                Tyiic           = repmat(Tyiic, length(ind), 1);
                Txiic           = repmat(Txiic, length(ind), 1);
                yyRc            = repmat(Ry(ind), 1, size(Tyiic, 2)) - Tyiic;
                xxRc            = repmat(Rx(ind), 1, size(Txiic, 2)) - Txiic;
                RRc             = sqrt(yyRc .^ 2 + xxRc .^ 2);
                
                ERc             = 1j * k / 4 * besselh(1, 1, -k * RRc) ./ (RRc + eps);
                
                ERcx            = -yyRc .* ERc;
                ERcy            = +xxRc .* ERc;
                
                ERtang          = ERcx .* cos(Rrad(ind) + pi / 2) + ERcy .* sin(Rrad(ind) + pi / 2);
                
                Enor            = FInc(calind, ii) / ERtang(calind);
                EPhi            = Enor * EPhi;
                Einc(:, ii)     = EPhi(:);
            end
            
            if para.disflag
                calind = 1 : length(FInc);
                ERtang = Enor * ERtang;
                figure;
                subplot(2, 2, 1.5)
                plot(calind, abs(FInc(calind,ii)), 'r-o');
                hold on;
                plot(calind, abs(ERtang(calind)), 'b--*');
                xlim([1 length(FInc)])
                subplot(2, 2, 3)
                plot(calind, real(FInc(calind,ii)), 'r-o');
                hold on;
                plot(calind, real(ERtang(calind)), 'b--*')
                xlim([1 length(FInc)])
                subplot(2, 2, 4)
                plot(calind, imag(FInc(calind,ii)), 'r-o');
                hold on;
                plot(calind, imag(ERtang(calind)), 'b--*')
                xlim([1 length(FInc)])
                plottools('on')
            end
            %             norm(ERtang(calind)-FInc(calind,ii))/norm(FInc(calind,ii))
            
        case PT.TM
            
            [NY, NX, ~]	    = meshgrid(ypuni(ny), xpuni(nx), zpuni(nz));
            
            Dx              = repmat(Rxuniq, 1, numel(NX)) - repmat(NX(:).', numel(Rxuniq), 1);
            Dy              = repmat(Ryuniq, 1, numel(NY)) - repmat(NY(:).', numel(Ryuniq), 1);
            R               = sqrt(Dx .^ 2 + Dy .^ 2);
            Phi	            = omega / 4 * besselh(0, 1, -k * R) * omega ^ 2 / (-1j * omega);
            
            % modified linear sampling method
            if isfield(para, 'nK') && para.nK > 0
                pars_LSM.PhiSin = cell(1,para.nK);
                pars_LSM.PhiCos    	= cell(1,para.nK);
                theta           = acos(Dx./R);
                for n = 1 : para.nK
                    thetan          = n * theta;% 13-07-2017
                    Phin            = omega/4*besselh(n,1,-k*R)*omega^2/(-1j*omega);
                    pars_LSM.PhiSin{n} 	= Phin.*(sin(thetan)); % 13-07-2017
                    pars_LSM.PhiCos{n} 	= Phin.*(cos(thetan)); % 13-07-2017
                end
            end
            
            k               = k / m_unit;
            Einc            = zeros(prodN, para.NTX);
            
            for ii = 1 : para.NTX
                Txiic           = Tx(NRX * ii - NRX + 1, :);
                Tyiic           = Ty(NRX * ii - NRX + 1, :);
                TyiiC           = repmat(Tyiic, size(NX));
                TxiiC           = repmat(Txiic, size(NX));
                yyRC            = NY * m_unit - TyiiC;
                xxRC            = NX * m_unit - TxiiC;
                RRC             = sqrt(yyRC .^ 2 + xxRC .^ 2);
                
                EPhi            = omega / 4 * besselh(0, 1, -k * RRC) * omega^2 / (-1j * omega);
                
                ind             = (1 : NRX) + (ii - 1) * NRX;
                Tyiic           = repmat(Tyiic, length(ind), 1);
                Txiic           = repmat(Txiic, length(ind), 1);
                yyRc            = repmat(Ry(ind), 1, size(Tyiic, 2)) - Tyiic;
                xxRc            = repmat(Rx(ind), 1, size(Txiic, 2)) - Txiic;
                RRc             = sqrt(yyRc .^ 2 + xxRc .^ 2);
                
                ER              = omega / 4 * besselh(0, 1, -k * RRc) * omega^2 / (-1j * omega);
                
                Enor            = sum(Finc(calind, ii)) / sum(ER(calind));
                EPhi            = Enor * EPhi;
                Einc(:, ii) 	= EPhi(:);
            end
            if para.disflag
                calind  = 1 : length(ER);
                ER      = Enor * ER;
                figure;
                subplot(2, 2, 1.5)
                plot(calind, abs(Finc(calind, ii)), 'r-o');
                hold on; plot(calind, abs(ER(calind)), 'b--*');
                xlim([1 length(ER)])
                subplot(2, 2, 3)
                plot(calind, real(Finc(calind, ii)), 'r-o');
                hold on;
                plot(calind, real(ER(calind)), 'b--*')
                xlim([1 length(ER)])
                subplot(2, 2, 4)
                plot(calind, imag(Finc(calind, ii)), 'r-o');
                hold on;
                plot(calind, imag(ER(calind)), 'b--*')
                xlim([1 length(ER)])
                plottools('on')
            end
            %             norm(ER(calind)-Finc(calind, ii))/norm(Finc(calind, ii))
    end
    Phi = full(Phi);
    %     save(Phiname, 'Phi', 'osc', 'grid3d', 'pars');
end



Phi             = Phi / osc.in_omega0();
Ainv         	= Ainv / osc.in_omega0() / osc.in_omega0();
% ii              = 1 : pars.Ninv(1);
% ii              = ii - pars.Ninv(1)/2;
% jj              = 1 : pars.Ninv(2);
% jj              = jj - pars.Ninv(2)/2;
% [II, JJ]        = meshgrid(ii, jj);
% RR              = sqrt(II.^2+JJ.^2);
% rr              = ceil(1.5*0.1/m_unit);
% Md              = RR > rr;














Md              = ones(Ninv(1), Ninv(2));

Md(6 : end - 5, 6 : end - 5) = 0;

Md              = Md(:);
rd              = InvInd(Md, para.pt);
pars.rd         = rd;

rchi            = InvInd(Md, PT.TM);
pars.rchi       = rchi;
N               = grid3dinv.N;
MdII            = reshape(Md, N(1), N(2));
vMdII           = MdII(:, round(N(2) / 2));
nxx             = find(diff(vMdII) ~= 0);
vMdII           = MdII(round(N(1) / 2), :);
nyy             = find(diff(vMdII) ~= 0);
MdII            = ones(size(MdII));

MdII(nxx(1) + 4 : nxx(2) - 4, nyy(1) + 4 : nyy(2) - 4) = 0;

rchiTV          = InvInd(MdII(:), PT.TM);
pars.rchiTV     = rchiTV;
pars.Nsrc       = para.NTX;
pars.Rindex     = Rindex;
pars.Runiq      = Runiq;

end