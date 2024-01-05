% time: 2023-5-2
% code by Baolei Wei

% noisy scenario:
% validating the numerical errors 
% on parameter estimation performance 

clc
clear
close all

% load routines
addpath('./utils')

% pars configuration
samplingManner = 0;                % 1: regular; 0: irregular
switch samplingManner
    case 1
        hOn = [0.35, 0.14, 0.07];  % interval
    otherwise
        hOn = [19, 49, 99];        % irregular
end

for K = 1:length(hOn)
    if samplingManner               % regularly sampling 
        h = hOn(K);
        ts = (h:h:7)';
    else                            % irregularly sampling 
        nn = hOn(K);
        rng(1);
        ts = 7*sort(unique([rand(nn,1);1] ) ); % nn = [19, 49, 99]
    end
    
    a1 = 1.6; a2 = -0.5; r = 1.5; y0 = 1.2;
    par.true = [a1 a2 r y0];
    
    % closed-form solution of ide model
    xideSol = @(t,p)(p(1)*p(4)^(1-p(3)) +p(2) ).*exp(p(1)*(1-p(3))*(t-t(1)) ).*...
                     ((p(4)^(1-p(3)) +p(2)/p(1) )*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(p(3)/(1-p(3) ) );
    xide = xideSol(ts,par.true); 
    
    yodeSol = @(t,p)((p(4)^(1-p(3) ) +p(2)/p(1) ).*exp(p(1)*(1-p(3)).*(t-t(1)) ) -p(2)/p(1) ).^(1/(1-p(3)) );

    %% measurement noise
    for nvr = [5 10 15]/100  % 5 10 15
        nrep = 500;
        pars = {'ode_sls', 'ode_nls', 'ide_sls', 'ide_nls'; ...
                nan(nrep,4), nan(nrep,4), nan(nrep,4), nan(nrep,4); ... % parameters 
                nan(nrep,1), nan(nrep,1), nan(nrep,1), nan(nrep,1)};    % evaluation metrics
      
        for irep = 1:nrep
            % data generation 
            rng(irep)
            xn = xide + nvr*std(xide)*randn(size(xide));
            
            % options for LM algorithm
            opts = optimoptions('lsqnonlin',...
                'Algorithm','levenberg-marquardt', ...
                'MaxIterations',8e6, 'MaxFunctionEvaluations',8e6);

            % initialisation of sls
            ryo = [2.0    0];           % initial
            ryl = [1.2 -inf];           % lower bound
            ryu = [3.0  inf];           % upper bound

            % initialisation of nls
            pl = [-inf -inf 1.2 -inf];  % lower bound
            pu = [ inf  inf 3.0  inf];  % upper bound

            %% ode-based model
            dt = diff(ts);              % physics-preserving Cusum operator
            yn = [0; cumsum(dt.*(xn(1:end-1)+xn(2:end)))/2];

            % sls for ODE
            lossfcn = @(ry)odeloss(ry,ts,yn + ry(2));
            ry.ode = lsqnonlin(lossfcn,ryo,ryl,ryu,opts);
            par.ode_sls = [a12' ry.ode];
               
            % nls for ODE
            po = par.ode_sls;           % initial guess for nls
            lossfcn = @(p) yodeSol(ts,p) - (yn+p(4));
            par.ode_nls = lsqnonlin(lossfcn,po,pl,pu,opts);

            %% ide-based model
            % sls for ide
            lossfcn = @(ry)ideloss(ry, ts, xn);
            ry.ide = lsqnonlin(lossfcn,ryo,ryl,ryu,opts);
            par.ide_sls = [a12' ry.ide];

            % nls for ide
            po = par.ide_sls;           % initial guess for nls
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
            mape.ide_sls = mean(abs(xf.ide_sls - xf.true)./xf.true)*100;
            mape.ode_sls = mean(abs(xf.ode_sls - xf.true)./xf.true)*100;
            mape.ide_nls = mean(abs(xf.ide_nls - xf.true)./xf.true)*100;
            mape.ode_nls = mean(abs(xf.ode_nls - xf.true)./xf.true)*100;            
            
            %% parmeters and metrics
            pars{2,1}(irep,:) = par.ode_sls; pars{2,2}(irep,:) = par.ode_nls;
            pars{2,3}(irep,:) = par.ide_sls; pars{2,4}(irep,:) = par.ide_nls;
            
            pars{3,1}(irep,:) = mape.ode_sls; pars{3,2}(irep,:) = mape.ode_nls;
            pars{3,3}(irep,:) = mape.ide_sls; pars{3,4}(irep,:) = mape.ide_nls;    
        end
        
%         %% significant test             
%         if samplingManner  % regular sampling
%             sigTestMatrix = [h*ones(500,1), nvr*ones(500,1), cell2mat(pars(3,1:4))];    
% 			writematrix(sigTestMatrix,'smape-regular.xlsx','WriteMode','append')
% 		else  % irregular sampling
%             sigTestMatrix = [nn*ones(500,1), nvr*ones(500,1), cell2mat(pars(3,1:4))];    
% 			writematrix(sigTestMatrix,'smape-irregular.xlsx','WriteMode','append')
% 		end
        
        %% tabs
        outputs = { samplingManner*ones(8,1), length(ts)*ones(8,1), nvr*ones(8,1), ...
                    [ mean(pars{2,1}) mean(pars{3,1}); std(pars{2,1}) std(pars{3,1}); ... % ode_sls
                      mean(pars{2,2}) mean(pars{3,2}); std(pars{2,2}) std(pars{3,2}); ... % ode_nls
                      mean(pars{2,3}) mean(pars{3,3}); std(pars{2,3}) std(pars{3,3}); ... % ide_sls
                      mean(pars{2,4}) mean(pars{3,4}); std(pars{2,4}) std(pars{3,4})] };  % ide_nls
        writematrix(cell2mat(outputs),'Noisy_pars.xlsx','WriteMode','append')
        
        %% figs
        f1 = figure;
        tiledlayout(1,5);
        clr = [255, 163,   0, 255;        % brown
               255,   0, 255, 255;        % cyan
                 0, 255, 255, 255;        % blue
                 0, 255,   0, 255]'/256;  % green 
             
        for i=1:4
            nexttile
            gboxplot({pars{2,1}(:,i), pars{2,2}(:,i), pars{2,3}(:,i), pars{2,4}(:,i)}, 1, clr)
            ytickformat('%.1f')
            switch i
                case 1
                    hold on; yline(a1,'-.b','linewidth',1); hold off
                    title('$\hat{a}_1$', 'interpreter','latex')
                case 2
                    hold on; yline(a2,'-.b','linewidth',1); hold off
                    title('$\hat{a}_2$', 'interpreter','latex')
                case 3
                    hold on; yline(r,'-.b','linewidth',1); hold off
                    title('$\hat{r}$', 'interpreter','latex')
                case 4
                    hold on; yline(y0,'-.b','linewidth',1); hold off
                    title('$\hat{\eta}_y$', 'interpreter','latex')
            end
            set(gca,'xticklabel', {''},'fontsize',13)
            grid on
        end

        nexttile
        gboxplot(pars(3,:), 1, clr)
        ytickformat('%.0f'); grid on 
        set(gca,'xticklabel', {''},'fontsize',13)
        title('$\mathrm{SMAPE}$','Interpreter','latex','fontsize',11)
        set(f1,'position',[200 300 850 180],'PaperSize',[15.5 5])  % 16 8
    
        if samplingManner
            exportgraphics(f1,'./figs/nvr_'+string(length(ts))+'-'+string(nvr*100)+'_regular.pdf')
        else
            exportgraphics(f1,'./figs/nvr_'+string(length(ts))+'-'+string(nvr*100)+'_irregular.pdf')
        end
    end    
end

