
clc
clear 
close all

%% data preprocessing 
R = readtable('trafficFlow.csv');
[y,m,d] = ymd(table2array(R(:,1)));
[h,~,~] = hms(table2array(R(:,1)));
S = [y,m,d,h, table2array(R(:,2:3))];    % junction vehicles

bt = 6;
for j = 1:4
    T{1, j} = S(S(:,5) == j,[1:4 6]);    % delete junction number
    J{1, j} = T{1, j}(bt:end,:);

    V{1, j} = reshape(J{1, j}(1:end-24-1+bt,1),24,[]);  % year
    V{2, j} = reshape(J{1, j}(1:end-24-1+bt,2),24,[]);  % month
    V{3, j} = reshape(J{1, j}(1:end-24-1+bt,3),24,[]);  % day
    V{4, j} = reshape(J{1, j}(1:end-24-1+bt,4),24,[]);  % hour
    V{5, j} = reshape(J{1, j}(1:end-24-1+bt,5),24,[]);  % vehicles
end

for j=1:4
    date.all = datetime(V{1, j},V{2, j},V{3, j},'Format','yyyy-MM-dd');
    date.num = weekday(date.all(1,:) );
    for wkd=1:7
        Date{wkd,j} = date.all(1, date.num == wkd);    % sunday - saturday
        Vehs{wkd,j} = V{5, j}(:, date.num == wkd);     % sunday - saturday
    end
end

junc_num = 1;
weeks = 50+[1:5];

% FUgure 
for wkd = 1:7
    figure(wkd)
    plot(Vehs{wkd, junc_num}(:,weeks),'LineWidth',2)
    hold on
    X(:,wkd) = mean(Vehs{wkd, junc_num}(:,weeks),2);
    plot(X(:,wkd),'LineWidth',2,'Color','k')
    xlim([0.5 24.5])
    grid on; hold off
    [~, dayname{wkd}] = weekday(Date{wkd,junc_num}(1,1),'long'); 
    display(string(Date{wkd,junc_num}(1,weeks(1))) +' ~ '+ string(Date{wkd,junc_num}(1,weeks(end))))
    title(string(dayname{wkd}))
    set(gca,'fontsize',22,...
        'XTick',[1 5 9 13 17 21 24],'XTickLabel',...
        {'5:00','9:00','13:00','17:00','21:00','1:00','4:00'})
    xtickangle(45)
    set(gcf,'position',[200 300 520 450]) 
%     exportgraphics(gcf,'figs/tflow-'+string(dayname{wkd})+'.pdf')
%     if wkd == 7
%         figure(8)
%         plot(Vehs{wkd, junc_num}(:,weeks),'LineWidth',2)
%         hold on
%         X(:,wkd) = mean(Vehs{wkd, junc_num}(:,weeks),2);
%         plot(X(:,wkd),'LineWidth',2,'Color','k')
%         xlim([0.5 24.5])
%         set(gca,'fontsize',22,...
%             'XTick',[1 5 9 13 17 21 24],'XTickLabel',...
%             {'5:00','9:00','13:00','17:00','21:00','1:00','4:00'})
%         xtickangle(45)
%         legend('1^{st} week','2^{nd} week',...
%             '3^{rd} week','4^{th} week','5^{th} week','average',...
%             'Orientation','vertical','Location','east','fontsize',22)
%         set(gcf,'position',[200 300 520 450]) 
%         exportgraphics(gcf,'tflow-legend.pdf')
%     end    
end

%% pattern dientification 
addpath('./utils')

ts = 1*(1:24)';

% closed-form solution of ide model
xideSol = @(t,p)(p(1)*p(4)^(1-p(3)) +p(2) ).*exp(p(1)*(1-p(3))*(t-t(1)) ).*...
                 ((p(4)^(1-p(3)) +p(2)/p(1) )*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(p(3)/(1-p(3) ) );
yodeSol = @(t,p)((p(4)^(1-p(3) ) +p(2)/p(1) ).*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(1/(1-p(3)) );


for wkd = 1:7
    xn = X(:, wkd);   

    % initialisation of sls
    ryo = [1.5    0];   % initial guess
    ryl = [1.2 -inf];   % lower bound
    ryu = [3.5  inf];   % upper bound
    
    % initialisation of nls
    pl = [-inf -inf 1.2 -inf];  % lower bound
    pu = [ inf  inf 3.5  inf];  % upper bound
    
    % options for LM algorithm
    opts = optimoptions('lsqnonlin', ...
        'Algorithm','levenberg-marquardt', ...
        'MaxIterations',8e9, 'MaxFunctionEvaluations',8e9);

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
    tf = ts;
    xf.ide_sls = xideSol(tf, par.ide_sls);
    xf.ide_nls = xideSol(tf, par.ide_nls);
    
    yf.ode_sls = yodeSol(tf, par.ode_sls);
    yf.ode_nls = yodeSol(tf, par.ode_nls);  
    xf.ode_sls = par.ode_sls(1)*yf.ode_sls + par.ode_sls(2)*yf.ode_sls.^par.ode_sls(3);
    xf.ode_nls = par.ode_nls(1)*yf.ode_nls + par.ode_nls(2)*yf.ode_nls.^par.ode_nls(3);

    %% metrics
    mape.ode_sls = mean(abs(xf.ode_sls - xn)./xn)*100;
    mape.ode_nls = mean(abs(xf.ode_nls - xn)./xn)*100;
    mape.ide_sls = mean(abs(xf.ide_sls - xn)./xn)*100;
    mape.ide_nls = mean(abs(xf.ide_nls - xn)./xn)*100;

    
    % Figure
    figure(10+wkd)
    plot(tf,xf.ode_sls,'-r', tf,xf.ode_nls,'-.g', ...
         tf,xf.ide_sls,':b', tf,xf.ide_nls,'--m', ...
         ts,xn,'.-k',  'markersize',25, 'LineWidth',2)
    xlim([0.5 24.5]); grid on
    title(string(dayname{wkd}))
    set(gca,'fontsize',22,...
        'XTick',[1 5 9 13 17 21 24],'XTickLabel',...
        {'5:00','9:00','13:00','17:00','21:00','1:00','4:00'})
    xtickangle(45)
    set(gcf,'position',[200 300 520 450]) 
%     exportgraphics(gcf,'figs/tfRES-'+string(dayname{wkd})+'.pdf')
%     if wkd == 7
%         figure(10+8)
%         plot(tf,xf.ode_sls,'-r', tf,xf.ode_nls,'-.g', ...
%              tf,xf.ide_sls,':b', tf,xf.ide_nls,'--m', ...
%              ts,xn,'.-k',  'markersize',25, 'LineWidth',2)
%         xlim([0.5 24.5]); grid on
%         set(gca,'fontsize',22,...
%             'XTick',[1 5 9 13 17 21 24],'XTickLabel',...
%             {'5:00','9:00','13:00','17:00','21:00','1:00','4:00'})
%         xtickangle(45)
%         legend('ODE_{sls}', 'ODE_{nls}', ...
%            '  IDE_{sls}', '  IDE_{nls}','average', ...
%            'Location','east', 'Orientation','vertical','fontsize',22)
%         xtickangle(45)
%         set(gcf,'position',[200 300 520 450]) 
%         exportgraphics(gcf,'tfRES-legend.pdf')
%     end    

    % outputs 
    pars_mape = [wkd par.ode_sls mape.ode_sls par.ide_sls mape.ide_sls; 
                 nan par.ode_nls mape.ode_nls par.ide_nls mape.ide_nls];
%     writematrix(pars_mape,'0-pars-mapes.xlsx','WriteMode','append')

    outputs= [xn xf.ode_sls abs(xf.ode_sls - xn)./xn*100 ...
                 xf.ode_nls abs(xf.ode_nls - xn)./xn*100 ...
                 xf.ide_sls abs(xf.ide_sls - xn)./xn*100 ...
                 xf.ide_nls abs(xf.ide_nls - xn)./xn*100; 
             nan nan mape.ode_sls nan mape.ode_nls ...
                 nan mape.ide_sls nan mape.ide_nls ];
%     writematrix(outputs,string(wkd)+'-'+string(dayname{wkd})+'.xlsx')
    
end





















