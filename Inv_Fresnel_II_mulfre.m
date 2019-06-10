clearvars; restoredefaultpath; add_path
%% load data
% TM
% str = 'FoamDielIntTM'; fre = 2 : 10; para.pt = PT.TM; NTX = 8; centre = [-5, 0, 0];
% str = 'FoamDielExtTM'; fre = 2 : 2 : 10; para.pt = PT.TM; NTX = 8; centre = [-20 0, 0];
% str = 'FoamTwinDielTM'; fre = 2 : 10; para.pt = PT.TM; NTX = 18; centre = [0 0, 0];
% str = 'FoamMetExtTM'; fre = 2 : 18; para.pt = PT.TM; NTX = 18; centre = [0 0, 0];
% TE
% str = 'FoamDielIntTE'; fre = 2 : 10; para.pt = PT.TE; NTX = 8; centre = [-5, 0, 0];
str = 'FoamDielExtTE'; fre = 2 : 2 : 10; para.pt = PT.TE; NTX = 8; centre = [-20 0, 0];
% str = 'FoamTwinDielTE'; fre = 2 : 10; para.pt = PT.TE; NTX = 18; centre = [0 0, 0];
% str = 'FoamMetExtTE'; fre = 2 : 18; para.pt = PT.TE; NTX = 18; centre = [0 0, 0];

xs              = 100; 
ys              = 100; 
para.m_unit     = 1.5e-3;
regSize         = [-3.0, 3.0; -3.0, 3.0]; % 12GHz
TxInterval      = 360 / NTX;
para.NTX        = NTX;
rho             = 0.20066;
DTr             = rho / 0.6964;
% DTr             = zeros(1, 18);

% if para.pt == PT.TE; DTr(2) = 0; end

para.Tr      	= 1.67 + DTr / 2;
para.Rr       	= 1.67 + DTr / 2;
rawdat          = load([str '.txt']);
rawdat(:, 1)    = TxInterval * (rawdat(:, 1) - 1) + 0.5;

%% Configuration

[Emea, Phi, A, chiinv, vJ, Einc, eTotInv, grid3d, pars] ...
    = Pre_InvFresnel_Conf(rawdat, fre, regSize, centre, xs, ys, para);

%% Display

run show_InitGuess.m

%% Iteration starts

rf                  = 1 : length(fre);      % select the frequency components for inversion
pars.bgchi          = 1;
pars.disflag        = true;                 % show the output every N iterations?
pars.terflag        = false;                %
pars.recflag        = false;                % save the records?
pars.omega          = pars.omega(rf);

pars.itenum       	= 16;                   % set iteration number

pars.recstr	= ['RecMF' str '_CCCSI'];

[mchi1, time1, chi1, vJ1, eTot1] = CCCSI_BrentMF(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

pars.recstr	= ['RecMF' str '_MRCSI'];

[mchi2, time2, chi2, vJ2, eTot2] = MRCSI_BrentMF(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

pars.recstr	= ['RecMF' str '_CSI'];

[mchi3, time3, chi3, vJ3, eTot3] = CSI_MF(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

%% save the inversion output

% save(['OMF' str  '.mat'], 'xinv', 'yinv', 'invdom', 'omega',...
%     'mchi1', 'time1', 'chi1', 'vJ1', 'eTot1', ...
%     'mchi2', 'time2', 'chi2', 'vJ2', 'eTot2', ...
%     'mchi3', 'time3', 'chi3', 'vJ3', 'eTot3', ...
%     'grid3d', 'ii', 'pars')

%% Show the inversion output

run show_output.m

