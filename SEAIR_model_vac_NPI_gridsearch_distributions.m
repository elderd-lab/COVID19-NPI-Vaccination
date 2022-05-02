% SEAIR MODEL - vac & NPI
%   SEAIR epidemic without births or deaths.
%   It contains an asymptomatic phase, A, vaccination, and a
%   portion of the population, alpha, 
%   practices non-pharmaceutical interventions (NPIs)
%
% Parameters drawn from distributions from published literature


% Non-varying parameters
b_V = 1; % relative transmission of asymptomatic
a = 0.35; % Asymptomatic parameter
xi = 0.95; % vaccination efficacy
epi = 0.50; % efficacy of NPI [Could vary as well - 30% to 70% per Ngonhala]
% 2% of pop vaccinated in Saad-Roy et al 2021 = 0.02
rho_avg = 0.002; % weighted average of vaccination rate
    
% Initial Conditions
N0 = 2.5e4;
Z0 = 0; % Initial Recovered
Z0Vac = 0; % Initial Vaccinated
Y0 = 25; % Initial Infected
V0 = 0; % Initial Asypmtomatic
W0 = 0; % Initial Exposed
X0 = N0 - (W0 + V0 + Y0 + Z0); % Initial # Susceptible

% Time Span
tend = 300; 
Tm = [0 tend];

% Parameters to be searched over
% alpha = 0.5; % proportion of pop practicing NPI
% rho = 0.002; % vaccination rate of individuals not practicing NPIs
% rho_NPI = 0.002; % vaccination rate of individuals practicing NPIs

alpha = 0.0:0.01:1.0; % proportion of pop practicing NPI
ratio = 0.10:0.05:5.0; % ratio of non_NPI/NPI vaccination rate

% Number of iterations
iter = 2e2;

% Storage matrices
[a1, b1] = size(alpha);
[a2, b2] = size(ratio);
Ymax = zeros(b1, b2, iter);
Zmax = zeros(b1, b2, iter);

for z = 1:iter % iteration loop
    % Default parameters w/ distribution draws
    %  NOTE: Used normal distribution to back calculate SIGMA. 
    %        95% = MU +/- 1.96*SIGMA.
    %        Log-normals upper confidence bounds were too wide.
    
    % mu = 2.7; 95% CI (1.6, 3.9) R0 - Davies et al - Lancet Public Health
    % Mean Runs - 2.7, 5.0, 10.0
    % Std Runs - 0.5612, 0.5612/1.5, 0.5612*1.5 
    R0 = normrnd(2.7, 0.5612*1.5); 
    % mu = 5.79; 95% CI (5.48, 6.11) Latency P - Brauner, Nature, Tab.S2
    gamma_W = 1/normrnd(5.79, 0.1582); 
    % mu = 3.45; 95% CI (3.24, 3.70) Infectious P - Li, Science
    gamma_V = 1/normrnd(3.45, 0.1071); 
    gamma_Y = 1/normrnd(3.45, 0.1071); % Same for asymptomatic and symptomatic
    beta = gamma_V*R0; % NOTE - needs to change if gamma_V != gamma_Y 

    for r = 1:b1 % alpha loop (NPI Compliance)

     % Initial conditions based on alpha
     X_initial=(1 - alpha(r)).*X0; W_initial=(1 - alpha(r)).*W0; 
     V_initial=(1 - alpha(r)).*V0; Y_initial=(1 - alpha(r)).*Y0;
     X_NPI_initial=alpha(r).*X0; W_NPI_initial=alpha(r).*W0; 
     V_NPI_initial=alpha(r).*V0; Y_NPI_initial=alpha(r).*Y0; 
     Z_initial=Z0; ZVac_initial = Z0Vac;    
    
      for c = 1:b2 % rho loop (Vaccination rate)

        % Calculate rho & rho_NPI based on ratio & weighted average
         rho_NPI = rho_avg/(alpha(r) + (1 - alpha(r)).*ratio(c));
         rho = ratio(c).*rho_NPI;
        
        % Call the ODE function 
        [t, pop]=ode45(@Diff_XWVYZ_vacNPI, Tm,[X_initial X_NPI_initial W_initial ...
        W_NPI_initial V_initial V_NPI_initial Y_initial Y_NPI_initial Z_initial ZVac_initial],...
        [],[beta epi b_V rho xi rho_NPI gamma_W a gamma_V gamma_Y]);
    
        X=pop(:,1); W=pop(:,3); V=pop(:,5); Y=pop(:,7);
        X_NPI=pop(:,2); W_NPI=pop(:,4); V_NPI=pop(:,6); Y_NPI=pop(:,8); 
        Z=pop(:,9); ZVac = pop(:,10);
        Ymax(r,c,z) = max(Y + Y_NPI); Zmax(r,c,z) = max(Z);
      end % rho (vaccine ratio) end 
    end % alpha (NPI compliance) loop end
end % iteration loop end

% Calculate the quantiles
quantYmax = quantile(Ymax, [0.1 .50 0.9], 3);
quantZmax = quantile(Zmax, [0.1 .50 0.9], 3);

figure

Levels = 0:0.001:0.12;

h(1) = subplot(2,3,1);
contourf(ratio,alpha,quantYmax(:,:,1)/N0,20)
title('Ymax - infecteds (2.5% - Lower bound)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Ymax as % pop.');

h(2) = subplot(2,3,2);
contourf(ratio,alpha,quantYmax(:,:,2)/N0,20)
title('Ymax - infecteds (50% - Median)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Ymax as % pop.');

h(3) = subplot(2,3,3);
contourf(ratio,alpha,quantYmax(:,:,3)/N0,20)
title('Ymax - infecteds (97.5% - Upper bound)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Ymax as % pop.');

h(4) = subplot(2,3,4);
contourf(ratio,alpha,quantZmax(:,:,1)/N0,20)
title('Zmax - recovered (2.5% - Lower bound)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Zmax as % pop.');

h(5) = subplot(2,3,5);
contourf(ratio,alpha,quantZmax(:,:,2)/N0,20)
title('Zmax - recovered (50% - Median)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Zmax as % pop.');

h(6) = subplot(2,3,6);
contourf(ratio,alpha,quantZmax(:,:,3)/N0,20)
title('Zmax - recovered (97.5% - Upper bound)')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Zmax as % pop.');

%% Calculates the differential rates used in the integration
function dPop=Diff_XWVYZ_vacNPI(t,pop, parameter)
 
 beta=parameter(1); epi=parameter(2); b_V = parameter(3); 
 rho = parameter(4); xi = parameter(5); rho_NPI = parameter(6);
 gamma_W = parameter(7); a = parameter(8); gamma_V = parameter(9); 
 gamma_Y = parameter(10);
 
 X=pop(1); X_NPI = pop(2);
 W=pop(3); W_NPI = pop(4);
 V=pop(5); V_NPI = pop(6);
 Y=pop(7); Y_NPI = pop(8); 
 Z=pop(9); Z_vac = pop(10);
 N = sum(pop);
 
 dPop=zeros(9,1);
 
 dPop(1)= -beta.*X.* ...
     (Y + (1 - epi).*Y_NPI + b_V.*V + (1 - epi).*b_V.*V_NPI)./N ...
     - rho.*xi.*X; % X
 dPop(2)= -beta.*(1 - epi).*X_NPI.* ...
     (Y + (1 - epi).*Y_NPI + b_V.*V + (1 - epi).*b_V.*V_NPI)./N ...
     - rho_NPI.*xi.*X_NPI; % X_NPI
 dPop(3) = beta.*X.* ...
     (Y + (1 - epi).*Y_NPI + b_V.*V + (1 - epi).*b_V.*V_NPI)./N -...
     gamma_W.*W; % W
 dPop(4) = beta.*(1-epi).*X_NPI.* ...
     (Y + (1 - epi).*Y_NPI + b_V.*V + (1 - epi).*b_V.*V_NPI)./N -...
     gamma_W.*W_NPI; % W_NPI
 dPop(5)= a.*gamma_W.*W - gamma_V.*V; % V
 dPop(6)= a.*gamma_W.*W_NPI - gamma_V.*V_NPI; % V_NPI
 dPop(7)= (1 - a).*gamma_W*W - gamma_Y.*Y; % Y 
 dPop(8)= (1 - a).*gamma_W*W_NPI - gamma_Y.*Y_NPI; % Y_NPI 
 dPop(9)= gamma_V.*V + gamma_V.*V_NPI + gamma_Y.*Y + ...
    gamma_Y.*Y_NPI; % Z Recovered 
 dPop(10) = rho.*xi.*X + rho_NPI.*xi.*X_NPI; % Z Vaccinated
end