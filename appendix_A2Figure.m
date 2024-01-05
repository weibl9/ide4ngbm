% trajectories of ide and ode
% ode generates y(t)
% ide generates y(t) and x(t)

clear; close

%%
h = .1;
t1 = h; tn = 7;
ts = (t1:h:tn)';

a1 = 1; a2 = -1; r = 2; eta = 0.2; % \eta_y = 0.2

%% analytic solution    
yodeSol = @(t) ((eta^(1-r ) +a2/a1 ).*exp(a1*(1-r).*(t-t1) ) -a2/a1 ).^(1/(1-r) );
yode = yodeSol(ts);

% ODE -> IDE
xodeSol = matlabFunction(diff(sym(yodeSol)));    % obtained from inverse integral operator    
xode = xodeSol(ts);    



xideSol = @(t) (a1*eta^(1-r) +a2 ).*exp(a1*(1-r)*(t-t1) ).*...
    ((eta^(1-r) +a2/a1 )*exp(a1*(1-r).*(t-t1) ) -a2/a1 ).^(r/(1-r));
xide = xodeSol(ts);

% IDE -> ODE
yideSol = matlabFunction(int(sym(xideSol)));     % obtained from integral operator 
yide = yideSol(ts); 


% FIGURE 
fsol = figure;
colororder({'r','b'})

yyaxis left
plot(ts,yode,'r.','markersize',16); hold on
plot(ts,yide,'-r','linewidth',  1); hold off
ylabel('$y(t)$','interpreter','latex')

yyaxis right
plot(ts,xode,'b.','markersize',16); hold on
plot(ts,xide,'-b','linewidth',  1); hold off
ylabel('$x(t)$','interpreter','latex')

xlabel('$t$','interpreter','latex')
legend('$y_\mathrm{ode}(t)$','$y_\mathrm{ide}(t)$',...
    '$x_\mathrm{ode}(t)$','$x_\mathrm{ide}(t)$',...
    'interpreter','latex','location','best','fontsize',15)
set(gca,'fontsize',12); 

set(fsol,'position',[200 300 600 450],'PaperSize',[15 10])
% print(fsol, 'trajectoryAnalytic','-dpdf')

%% numerical solutions
dt = diff(ts);

% IDE -> ODE
yide_num = eta +[0; dt.*cumsum(xide(1:end-1)+xide(2:end))/2 ];

% ODE -> IDE 
xode_num = a1*yode + a2*yode.^r;

% FIGURE 
fnum = figure;
colororder({'r','b'})

yyaxis left
plot(ts,yode,'r.','markersize',16); hold on
plot(ts,yide_num,'-r','linewidth',  1); hold off
ylabel('$y(t)$','interpreter','latex')

yyaxis right
plot(ts,xode_num,'b.','markersize',16); hold on
plot(ts,xide,'-b','linewidth',  1); hold off
ylabel('$x(t)$','interpreter','latex')
set(gca,'fontsize',12)

xlabel('$t$','interpreter','latex')
legend('$y_\mathrm{ode}(t)$','$y_\mathrm{ide}(t)$',...
    '$x_\mathrm{ode}(t)$','$x_\mathrm{ide}(t)$',...
    'interpreter','latex','location','best','fontsize',15)

set(fnum,'position',[200 300 600 450],'PaperSize',[15 10])
% print(fnum, 'trajectoryNumeric','-dpdf')


% FIGURE
ferr = figure;
subplot(2,1,1)
plot(ts, yode-yide_num,'.-r','markersize',16)
ylabel('$$y_{\mathrm{ode}}(t)-y_{\mathrm{ide}}(t)$$','interpreter','latex')
set(gca,'fontsize',12);

subplot(2,1,2)
plot(ts, xode_num-xide,'.-b','markersize',16)
ylabel('$$x_{\mathrm{ode}}(t)-x_{\mathrm{ide}}(t)$$','interpreter','latex') 
set(gca,'fontsize',12)

xlabel('$t$','interpreter','latex')

set(ferr,'position',[200 300 600 450],'PaperSize',[15 10])
% print(ferr,'trajectoryNumericError','-dpdf')







