close all
clear
clc

%% implementing ode45 to describe Cg, Cc, Cp, and Ce
% Recall that Ce is just excretory and is included to provide the
% mass balance for summing across all regions.
F = 0.14; %CHANGE THIS FOR EACH DRUG, use for oral admin only
V_gi = 0.120; % use this for oral admin only
D = 20; %CHANGE THIS FOR EACH DRUG
V_c = 1960;

V_p = 3900; %CHANG THIS FOR EACH DRUG
Q = 1960 / 60; %(mL/hour) %CHANGE THIS FOR EACH DRUG
%Sometimes they use intercompartmental clearance Q, kc = q/vc, kp = q/vp
kc = Q/V_c;

kp = Q/V_p;


%now calculating absorbance and elimination
t_el = 14 * 60; %CHANGE THIS FOR EACH DRUG, half life
ke = log(2) / t_el;
ke;

%t_abs = 11.5; %CHANGE THIS FOR EACH DRUG, half life, use for oral admin
%ka = log(2) / t_abs; %use for oral admin

ka = 2.5 / 60; %use for iv admin or known ka
ka;

%ke_1 = ke * 1.3;
ke_1 = ke;
%ka_1 = ka * 0.6;
ka_1 = ka;

%C_g_i = D*F / V_gi;
C_g_i = 0;
C_iv = D/V_c;
%C_iv = 0;

x=[C_g_i C_iv 0 0]; %CHANGE THIS IF IV INJECTION
tspan = 0:4000;
% 
% [t,y]=ode45(@odefun,tspan,x);
% P = y(:, 3);
% C = y(:, 2);


% Analytical Solution for one compartment IV
% By solving the DiffEQ we get the form C = C_0 * e^(-k_e * t)
% 
% C_0 = D/ (V_c + V_p); %see above
% t_as_IV = [0:4000];
% C_as_IV = [];
% y = 0;
% x = 0;
% 
% for i = 0:4000
%     y = C_0 * exp(-ke * x);
%     C_as_IV = [C_as_IV, y];
%     x = x + 1;
% end
% C_as_IV_Inj = C_as_IV;
% Avg_C_as_IV_Inj = mean(C_as_IV());
% 
% tspan = [0 4000];
% [t_ml_IV,C_ml_IV] = ode45(@(t_ml_IV,C_ml_IV) -ke * C_ml_IV, tspan, C_0);
% 
% 
% % % Analytical Solution for one compartment oral
% % % Solving the system of Differential Equations we get 
% % % C = (D * F * k_in) / (V_c * (k_in - k_e)) * (exp(-k_e * t) - exp(-k_in *
% % % t))
% % 
% % C_0 = 0; %see above
% % t_as_GI = [0:4000];
% % C_as_GI = [];
% % y = 0;
% % x = 0;
% % V_c1 = V_p + V_c;
% % %V_c1 = 0.9 * 70;
% % 
% % for i = 0:4000
% %     y = ((D * F * ka_1) / (V_c1 * (ka_1 - ke_1))) * (exp(-ke_1 * x) - exp(-ka_1 * x));
% %     C_as_GI = [C_as_GI, y];
% %     x = x + 1;
% % end
% % C_as_Oral = C_as_GI;
% % Avg_C_as_Oral = mean(C_as_Oral);
% % 
% % %matlab solving
% % C_0 = 0;
% % [t_ml_o,C_ml_o] = ode45(@(t_ml,C_ml) (((ka_1 * F * D * exp(-ka_1 * t_ml))/V_c1) - ((ke_1) * C_ml)), tspan, C_0);

% figure('name','ODE 1')
% %t_o = t + 5;
% plot(t, C)
% hold on
% plot(t, P)
% hold on
% %plot(t, C_as_Oral)
% plot(t_as_IV, C_as_IV)
% hold on
% %plot(t, C_ml_o)
% plot(t_ml_IV, C_ml_IV)
% 
% hold off
% xlabel('Time (minutes)')
% ylabel('Concentration (mg/L)')
% title('Time-dependent Concentration')
% axis('tight')
% legend('C-2C','P-2C', '1-C-as', '1-C-ml')

%% Introducing k_e values
% Parameters
k_values = [0.0167, 0.0833, 0.167];  % Different elimination rate constants (1/sec)
C_A = 0.1;          % Initial concentration of drug entering the compartment (mg/mL)
tspan = [0 200];    % Time span for simulation (seconds)

% Create a figure for plotting
figure;

% Initialize cell array to store legend labels
legend_labels = cell(length(k_values), 1);

% Loop through different values of k and plot results on the same graph
hold on; % Allow multiple plots on the same graph
for i = 1:length(k_values)
    k = k_values(i);
    
    % Solve the differential equation
    dCdt = @(t, C) -k * C;
    [t, C] = ode45(dCdt, tspan, C_A);
    
    % Plot concentration vs. time and store legend label
    plot(t, C, 'LineWidth', 2);
    legend_labels{i} = ['k = ', num2str(k), ' 1/sec'];
end

% Customize plot appearance
title('Drug Elimination with Varying Elimination Rates');
xlabel('Time (seconds)');
ylabel('Concentration (mg/mL)');
grid on;

% Add legend
legend(legend_labels);

hold off; % Stop allowing multiple plots on the same graph

%% Now extending to clearance rate values
Cl_values = [0.01, 0.05, 0.1];  % Different clearance rates (mL/sec)
C_A = 0.1;          % Initial concentration of drug entering the compartment (mg/mL)
tspan = [0 600];    % Time span for simulation (seconds)

% Create a figure for plotting
figure;

% Initialize cell array to store legend labels
legend_labels = cell(length(Cl_values), 1);

% Loop through different values of Cl and plot results on the same graph
hold on; % Allow multiple plots on the same graph
for i = 1:length(Cl_values)
    Cl = Cl_values(i);
    
    % Calculate k based on Cl and C_A
    k = Cl * C_A;
    
    % Solve the differential equation
    dCdt = @(t, C) -k * C;
    [t, C] = ode45(dCdt, tspan, C_A);
    
    % Plot concentration vs. time and store legend label
    plot(t, C, 'LineWidth', 2);
    legend_labels{i} = ['Cl = ', num2str(Cl), ' mL/sec'];
end

% Customize plot appearance
title('Renal Elimination with Varying Clearance Rates');
xlabel('Time (seconds)');
ylabel('Concentration (mg/mL)');
grid on;

% Add legend
legend(legend_labels);

hold off; % Stop allowing multiple plots on the same graph
%% Now adding extraction ratios
% Parameters
E_values = [0.1, 0.5, 0.9];  % Different extraction ratios (unitless)
Q = 16.67;          % Volumetric flow rate (mL/sec)
C_A = 0.1;          % Initial concentration of drug entering the compartment (mg/mL)
tspan = [0 10];    % Time span for simulation (seconds)

% Create a figure for plotting
figure;

% Initialize cell array to store legend labels
legend_labels = cell(length(E_values), 1);

% Loop through different values of E and plot results on the same graph
hold on; % Allow multiple plots on the same graph
for i = 1:length(E_values)
    E = E_values(i);
    
    % Calculate Cl based on E and Q
    Cl = E * Q;
    
    % Calculate k based on Cl and C_A
    k = Cl * C_A;
    
    % Solve the differential equation
    dCdt = @(t, C) -k * C;
    [t, C] = ode45(dCdt, tspan, C_A);
    
    % Plot concentration vs. time and store legend label
    plot(t, C, 'LineWidth', 2);
    legend_labels{i} = ['E = ', num2str(E)];
end

% Customize plot appearance
title('Drug Elimination with Varying Extraction Ratios (E)');
xlabel('Time (seconds)');
ylabel('Concentration (mg/mL)');
grid on;

% Add legend
legend(legend_labels);

hold off; % Stop allowing multiple plots on the same graph
