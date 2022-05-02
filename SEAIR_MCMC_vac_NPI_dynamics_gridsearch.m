% SEAIR MODEL with vaccination and NPI
%   SEAIR epidemic without births or deaths. 
%   Solves the deterministic skeleton (currently based on no. of individs.)
%   Runs model of stochastic dynamics (based on # of individuals)
%   X,W (Exposed),V (Asymptomatic), Y, and Z follows Keeling and Rohani 
%   book where tracking individuals.

% Default parameters w/ distribution draws
R0 = 2.7; % 95% CI (1.6, 3.9) Davies et al - Lancet Public Health;
gamma_W = 1/5.79; % 95% CI (5.48, 6.11) Latency P - Brauner, Nature, Tab.S2
gamma_V = 1/3.45; % 95% CI (3.24, 3.70) Infectious P - Li, Science
gamma_Y = 1/3.45; % Same for asymptomatic and symptomatic
beta = gamma_V*R0; % NOTE - needs to change if gamma_V != gamma_Y

% Non-varying parameters
b_V = 1; % relative transmission of asymptomatic
% [Note Sah et al estimate A as 35.1 (30.7, 39.9)]
a = 0.35; % 0.1 - Asymptomatic parameter 
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
tend = 300; 
Tm = [0 tend];

% Number of MCMC iterations
niter = 100;

% Parameters to be searched over
% alpha = 0.5; % proportion of pop practicing NPI
% rho = 0.002; % vaccination rate of individuals not practicing NPIs
% rho_NPI = 0.002; % vaccination rate of individuals practicing NPIs

alpha = 0.0:0.05:1.0; % 0.1 increment proportion of pop practicing NPI
ratio = 0.01:.25:5.0; % 0.05 increment ratio of non_NPI/NPI vaccination rate

[a1, b1] = size(alpha);
[a2, b2] = size(ratio);

% Worst-case scenario
%  NOTE: Mean of stochastic simulations will be the ODEs
Ymaxmax = zeros(b1,b2); % Max of max Y for stochastic simulation
Zmaxmax = zeros(b1,b2); 

for r = 1:b1 % alpha loop (NPI compliance)

    % Initial conditions based on alpha                      
    X_initial=(1 - alpha(r)).*X0; W_initial=(1 - alpha(r)).*W0; 
    V_initial=(1 - alpha(r)).*V0; Y_initial=(1 - alpha(r)).*Y0;
    X_NPI_initial=alpha(r).*X0; W_NPI_initial=alpha(r).*W0; 
    V_NPI_initial=alpha(r).*V0; Y_NPI_initial=alpha(r).*Y0; 
    Z_initial=Z0; ZVac_initial=Z0Vac;

    for c = 1:b2 % rho loop (Vaccination rate)

        % Calculate rho & rho_NPI based on ratio & weighted average
        rho_NPI = rho_avg/(alpha(r) + (1 - alpha(r)).*ratio(c));
        rho = ratio(c).*rho_NPI;
    
        % The main iteration 
        [t, pop]=ode45(@Diff_XWVYZ_vacNPI, Tm,[X_initial X_NPI_initial ...
        W_initial W_NPI_initial V_initial V_NPI_initial Y_initial ... 
        Y_NPI_initial Z_initial ZVac_initial],...
        [],[beta epi b_V rho xi rho_NPI gamma_W a gamma_V gamma_Y]);

        % The Stochastic Model - start array
        zstate = zeros(1, 11, 1); % time, [tm states], iteration

        % Transition Matrix
        % States:  X X_NPI W W_NPI V V_NPI Y Y_NPI Z
                                        % Transition Matrix
        lambda = [-1 0 1 0 0 0 0 0 0 0; % BXY/N
                  -1 0 1 0 0 0 0 0 0 0; % (1 - epi)*BXY_NPI/N
                  -1 0 1 0 0 0 0 0 0 0; % b_V*BXV/N
                  -1 0 1 0 0 0 0 0 0 0; % (1 - epi)*b_V*BXV_NPI/N
                  -1 0 0 0 0 0 0 0 0 1; % rho*xi*X
                   0 -1 0 1 0 0 0 0 0 0; % (1 - epi)*BX_NPI*Y/N
                   0 -1 0 1 0 0 0 0 0 0; % (1 - epi)^2*BX_NPI*Y_NPI/N
                   0 -1 0 1 0 0 0 0 0 0; % (1 - epi)*BX_NPI*V/N
                   0 -1 0 1 0 0 0 0 0 0; % (1 - epi)^2*BX_NPI*V_NPI/N
                   0 -1 0 0 0 0 0 0 0 1; % rho_NPI*xi*X_NPI
                   0 0 -1 0 1 0 0 0 0 0; % a*gamma_W*W
                   0 0 0 -1 0 1 0 0 0 0; % a*gamma_W*W_NPI
                   0 0 -1 0 0 0 1 0 0 0; % (1-a)*gamma_W*W
                   0 0 0 -1 0 0 0 1 0 0; % (1-a)*gamma_W*W_NPI
                   0 0 0 0 -1 0 0 0 1 0; % gamma_V*V
                   0 0 0 0 0 -1 0 0 1 0; % gamma_V*V_NPI
                   0 0 0 0 0 0 -1 0 1 0; % gamma_Y*Y
                   0 0 0 0 0 0 0 -1 1 0];% gamma_Y*Y_NPI
          
        % For loop for iterations
        for i = 1:niter
            tm = 0;
            tm_index = 1;
            vstate = [X_initial X_NPI_initial W_initial ...
                W_NPI_initial V_initial V_NPI_initial Y_initial ...
                Y_NPI_initial Z_initial ZVac_initial];

            % While loop for stochastic realizations
            while vstate(3) + vstate(4) + vstate(5) + vstate(6) ...
                    + vstate(7) + vstate(8) > 0 && tm <= tend
                zstate(tm_index, :, i) = [tm vstate];
                X = vstate(1);
                X_NPI = vstate(2);
                W = vstate(3);
                W_NPI = vstate(4);
                V = vstate(5);
                V_NPI = vstate(6);
                Y = vstate(7);
                Y_NPI = vstate(8);
                Z = vstate(9);
                ZVac = vstate(10);
                N = X + X_NPI + W + W_NPI + V + V_NPI + Y + Y_NPI + Z + ZVac;
                vec_p = [beta.*X.*Y./N (1-epi).*beta.*X.*Y_NPI./N ...
                b_V.*beta.*X.*V./N (1-epi).*b_V.*beta.*X.*V_NPI./N ...
                rho.*xi.*X ...
                (1-epi).*beta.*X_NPI.*Y./N (1-epi).^2*beta.*X_NPI.*Y_NPI./N ...
                (1-epi).*beta.*X_NPI.*V./N (1-epi).^2.*beta.*X_NPI.*V_NPI./N ...
                rho_NPI.*xi.*X_NPI a.*gamma_W.*W a.*gamma_W.*W_NPI ...
                (1 - a).*gamma_W.*W (1 - a).*gamma_W.*W_NPI ...
                gamma_V.*V gamma_V.*V_NPI gamma_Y.*Y gamma_Y.*Y_NPI];
            
                delta_t = 1/sum(vec_p);
                vec_l = poissrnd(vec_p.*delta_t, size(vec_p));
        
                vstate = vstate + vec_l*lambda;

                vstate(vstate < 0) = 0;
                tm = tm + delta_t;
                tm_index = tm_index + 1;
        
                %   YY = [round(sum(vstate))];
                %   disp(YY)
            end % Stochastic while loop end
        ZZ = ['Doing SDE realization: ', num2str(i),' of ', num2str(niter), ...
            ' vstate ' num2str(round(vstate)), ' total:', ...
            num2str(round(sum(vstate)))];
        disp(ZZ)
        end % Stochastic for loop niter end

% Reformating the Stochastic Results
[aa,bb,cc] = size(zstate);

mcmc_time = zstate(:,1,:);
mcmc_time = reshape(mcmc_time, aa, cc);

X_rep = zstate(:,2,:);
X_rep = reshape(X_rep, aa, cc);

X_NPI_rep = zstate(:,3,:);
X_NPI_rep = reshape(X_NPI_rep, aa, cc);

W_rep = zstate(:,4,:);
W_rep = reshape(W_rep, aa, cc);

W_NPI_rep = zstate(:,5,:);
W_NPI_rep = reshape(W_NPI_rep, aa, cc);

V_rep = zstate(:,6,:);
V_rep = reshape(V_rep, aa, cc);

V_NPI_rep = zstate(:,7,:);
V_NPI_rep = reshape(V_NPI_rep, aa, cc);

Y_rep = zstate(:,8,:);
Y_rep = reshape(Y_rep, aa, cc);

Y_NPI_rep = zstate(:,9,:);
Y_NPI_rep = reshape(Y_NPI_rep, aa, cc);

Z_rep = zstate(:,10,:);
Z_rep = reshape(Z_rep, aa, cc);

ZVac_rep = zstate(:,11,:);
ZVac_rep = reshape(ZVac_rep, aa, cc);

% Formatting  Deterministic Results
Xpop = pop(:,1); X_NPIpop = pop(:,2); Wpop = pop(:,3); W_NPIpop = pop(:,4); 
Vpop = pop(:,5); V_NPIpop = pop(:,6); Ypop = pop(:,7); Y_NPIpop = pop(:,8);
Zpop = pop(:,9); Zvac = pop(:,10);

% Store Ymaxmax Zmaxmax
Yrepmax = max(Y_rep + Y_NPI_rep, [], 2);
Ymaxmax(r,c) = max(Yrepmax);
Zrepmax = max(Z_rep + ZVac_rep, [], 2);
Zmaxmax(r,c) = max(Zrepmax);

ZZZ = ['Finished loop: ratio (vaccination rate), ', num2str(ratio(c)), ...
            ' alpha (NPI compliance), ', num2str(alpha(r))];
disp(ZZZ)

    end % for loop vaccination (rho) end
    
    if r == b1
        save NPIgridsearch.mat
    end
    
end % for loop NPI compliance (alpha) end

RunTime = toc;
% plots
figure
contourf(ratio,alpha,Ymaxmax/N0,20)
title('Stochastic Max of Ymax - infecteds')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Ymax as % pop.');

figure
contourf(ratio,alpha,Zmaxmax/N0,20)
title('Stochastic Max of Zmax - recovered')
ylabel('alpha, prop practicing NPI'); xlabel('ratio of Non-NPI/NPI Vac');
colormap('parula');
hcb = colorbar;
title(hcb,'Zmax as % pop.');


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
    gamma_Y.*Y_NPI; % Z
dPop(10) = rho.*xi.*X + rho_NPI.*xi.*X_NPI; % Z Vac
end


