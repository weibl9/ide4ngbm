% time: 2023-5-2
% code by Baolei Wei

% noise-free scenario:
%   validating the numerical errors 
%   on parameter estimation performance 

clc
clear
close all

%% load routines
addpath('./utils')

samplingManner = 1; 

switch samplingManner
    case 1       % regularly sampling 
        h = 0.7; 
        ts = (h:h:7)';
    otherwise    % irregularly sampling
        m = 9; rng(7)
        ts = 7*sort(unique([rand(m,1); 1] ) );
end

% true paramters
a1 = 1.6; a2 = -0.5; r = 1.5; y0 = 1.2;
par.true = [a1 a2 r y0];

% closed-form solution of ide model
xideSol = @(t,p)(p(1)*p(4)^(1-p(3)) +p(2) ).*exp(p(1)*(1-p(3))*(t-t(1)) ).*...
                 ((p(4)^(1-p(3)) +p(2)/p(1) )*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(p(3)/(1-p(3) ) );
xide = xideSol(ts,par.true); 

yodeSol = @(t,p)((p(4)^(1-p(3) ) +p(2)/p(1) ).*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(1/(1-p(3)) );

%% data generation and initialization
xn = xide;

% initialisation of sls
ryo = [2.0    0];   % initial guess
ryl = [1.2 -inf];   % lower bound
ryu = [3.0  inf];   % upper bound

% initialisation of nls
pl = [-inf -inf 1.2 -inf];  % lower bound
pu = [ inf  inf 3.0  inf];  % upper bound

% options for LM algorithm
opts = optimoptions('lsqnonlin', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxIterations',8e6, 'MaxFunctionEvaluations',8e6);


%% ode-based model    
dt = diff(ts);                  % physics-preserving Cusum operator
yn = [0; cumsum(dt.*(xn(1:end-1)+xn(2:end)))/2];

% sls for ODE
lossfcn = @(ry)odeloss(ry,ts,yn + ry(2));
ry.ode = lsqnonlin(lossfcn,ryo,ryl,ryu,opts);
par.ode_sls = [a12' ry.ode];
   
% nls for ODE
po = par.ode_sls;               % initial guess for nls
lossfcn = @(p) yodeSol(ts,p) - (yn+p(4));
par.ode_nls = lsqnonlin(lossfcn,po,pl,pu,opts);

[norm(lossfcn(po),2)^2 norm(lossfcn(par.ode_nls),2)^2]   % why nls<sls

%% ide-based model
% sls for ide
lossfcn = @(ry)ideloss(ry, ts, xn);
ry.ide = lsqnonlin(lossfcn,ryo,ryl,ryu,opts);
par.ide_sls = [a12' ry.ide];

% nls for ide
po = par.ide_sls;               % initial guess for nls
lossfcn = @(p)xideSol(ts, p) - xn;
par.ide_nls = lsqnonlin(lossfcn,po,pl,pu,opts);


%% trajectories  
tf = linspace(min(ts), max(ts), 7000)';
xf.true = xideSol(tf, par.true);
xf.ide_sls = xideSol(tf, par.ide_sls);
xf.ide_nls = xideSol(tf, par.ide_nls);

yf.ode_sls = yodeSol(tf, par.ode_sls);
yf.ode_nls = yodeSol(tf, par.ode_nls);  
xf.ode_sls = par.ode_sls(1)*yf.ode_sls + par.ode_sls(2)*yf.ode_sls.^par.ode_sls(3);
xf.ode_nls = par.ode_nls(1)*yf.ode_nls + par.ode_nls(2)*yf.ode_nls.^par.ode_nls(3);

%% metrics
mape.ode_sls = mean(abs(xf.ode_sls - xf.true)./xf.true)*100;
mape.ode_nls = mean(abs(xf.ode_nls - xf.true)./xf.true)*100;
mape.ide_sls = mean(abs(xf.ide_sls - xf.true)./xf.true)*100;
mape.ide_nls = mean(abs(xf.ide_nls - xf.true)./xf.true)*100;


%%  figure
mrk = 25;
lwd = 1.5;

f1 = figure;
plot(ts, xide, '.k','markersize',mrk); hold on
plot(tf, xf.ode_sls, '-r', 'linewidth', lwd)
plot(tf, xf.ode_nls, '--b', 'linewidth',lwd)
plot(tf, xf.ide_sls, ':g', 'linewidth', lwd)
plot(tf, xf.ide_nls, '-.k', 'linewidth',lwd); hold off
xlim([0 7.2]); grid on
xlabel('$t$','interpreter','latex')
ylabel('$x(t)$','interpreter','latex')
set(gca,'fontsize',15,'TickLabelInterpreter','latex')
legend('True',...
    'ODE$_\mathrm{sls}$','ODE$_\mathrm{nls}$',...
    '~IDE$_\mathrm{sls}$','~IDE$_\mathrm{nls}$',...
    'location','northeast', 'fontsize',13, 'interpreter','latex')
set(f1,'position',[10 10 650 400],'PaperSize',[15 10])


%% save figure
switch samplingManner
    case 1       % regularly sampling 
        exportgraphics(f1,'./figs/Noisefree_regular.pdf')
        fprintf('\n Regularly sampling: \n'); disp(par);
    otherwise    % irregularly sampling
        exportgraphics(f1,'./figs/Noisefree_irregular.pdf')
        fprintf('\n Irregularly sampling: \n'); disp(par);
end






