% SEAIR MODEL - vac & NPI
%   SEAIR epidemic without births or deaths.
%   (S)usceptible-(E)xposed-(A)symptomatic-(I)nfected-(R)ecovered
%   It contains an asymptomatic phase, A, vaccination, and a
%   portion of the population, alpha, 
%   practices non-pharmaceutical interventions (NPIs).
%

% Default parameters w/ distribution draws
R0 = 2.7; % 95% CI (1.6, 3.9) Davies et al - Lancet Public Health;
gamma_W = 1/5.79; % 95% CI (5.48, 6.11) Latency P - Brauner, Nature, Tab.S2
gamma_V = 1/3.45; % 95% CI (3.24, 3.70) Infectious P - Li, Science
gamma_Y = 1/3.45; % Same for asymptomatic and symptomatic
beta = gamma_V*R0; % NOTE - needs to change if gamma_V != gamma_Y

% Non-varying parameters
b_V = 1; % relative transmission of asymptomatic
a = 0.35; % Asymptomatic parameter
xi = 0.95; % vaccination efficacy
epi = 0.50; % efficacy of NPI
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
tend = 250; 
Tm = [0 tend];

% Parameters to be searched over
% alpha = 0.5; % proportion of pop practicing NPI
% rho = 0.002; % vaccination rate of individuals not practicing NPIs
% rho_NPI = 0.002; % vaccination rate of individuals practicing NPIs

alpha = [0.1 0.25 0.5]; % 0.0:0.01:1.0; % proportion of pop practicing NPI
ratio = 1; % 0.01:0.05:5.0; % ratio of non_NPI/NPI vaccination rate

[a1, b1] = size(alpha);
[a2, b2] = size(ratio);
Ymax = zeros(b1,b2);
Zmax = zeros(b1,b2);

% MODIFY DEPENDING ON STORAGE NEED or COMMENT OUT
PopStorage = -1*ones(220,11,3); % Time Series, time + Pop State, Compliance

figure; 

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
        Z=pop(:,9);
        Ymax(r,c) = max(Y + Y_NPI); Zmax(r,c) = max(Z);
    
    end
    
    [a1, b1] = size(pop);
    PopStorage(1:a1,:,r) = [t pop];
    
    subplot(1,3,r)
    h=plot(t,(W + W_NPI)/N0,'-g',t,(V + Y)/N0,'-k',t,(V_NPI + Y_NPI)/N0, '-r'); 
    xlabel 'Time';
    ylabel 'Exposed, Non-compliant, and Compliant';
    
end

% Calculates the differential rates used in the integration
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
 
 dPop=zeros(10,1);
 
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