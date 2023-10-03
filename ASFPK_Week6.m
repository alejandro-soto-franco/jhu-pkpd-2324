% Define the fixed excretion rate constant k_ex
k_ex = 0.01; % You can adjust this value as needed

% Define a range of metabolic rate constants k_m to study
k_m_values = [0.1, 0.5, 1.0]; % You can adjust these values as needed

% Time span for the simulation
tspan = [0, 24]; % Time span from 0 to 24 hours

% Create a logarithmic time span for the second plot
log_tspan = logspace(-2, log10(24), 1000);

% Create a figure for the concentration over time plot
figure;
subplot(2, 1, 1);
hold on;

% Loop through different k_m values and plot concentration over time for each
for i = 1:length(k_m_values)
    k_m = k_m_values(i);
    
    % Define the elimination rate constant ke as km + k_ex
    k_e = k_m + k_ex;

    % Define the ODE function
    dCdt = @(t, C) -k_e * C;

    % Initial concentration of dextromethorphan
    C0 = 100; % Initial concentration (you can adjust this as needed)

    % Solve the ODE using ode45
    [t, C] = ode45(dCdt, tspan, C0);

    % Plot concentration over time
    plot(t, C, 'DisplayName', sprintf('k_m = %.2f', k_m));
end

hold off;

xlabel('Time (hours)');
ylabel('Concentration (ng/mL)');
title('Drug Concentration over Time during under Metabolism');
legend('Location', 'Best');

% Create a logarithmic time axis plot
subplot(2, 1, 2);
hold on;

% Loop through different k_m values and plot log(concentration) versus time for each
for i = 1:length(k_m_values)
    k_m = k_m_values(i);
    
    % Define the elimination rate constant ke as km + k_ex
    k_e = k_m + k_ex;

    % Define the ODE function
    dCdt = @(t, C) -k_e * C;

    % Initial concentration of dextromethorphan
    C0 = 100; % Initial concentration (you can adjust this as needed)

    % Solve the ODE using ode45 with logarithmic time span
    [t, C] = ode45(dCdt, log_tspan, C0);

    % Plot log(concentration) versus time
    plot(t, log(C), 'DisplayName', sprintf('k_m = %.2f', k_m));
end

hold off;

xlabel('Time (hours)');
ylabel('log(Concentration)');
title('Logarithm of Drug Concentration over Time during under Metabolism');
legend('Location', 'Best');

%% Isoniazid Michaelis-Menten model
% Constants
Dose = 250; % mg
V_max = 30.87; % mg/mL/hr
K_m = 1.669; % mg/mL
R = 27.75; % hours
Saturation = 0.6; % 60% saturation

% Time span for simulation (e.g., 0 to 100 hours)
tspan = [0 100];

% Define a range of substrate concentrations
C_s_range = linspace(0, 100, 1000); % ng/mL

% Initialize arrays to store results
P_values = zeros(length(C_s_range), 1);
ReactionRate_values = zeros(length(C_s_range), 1);

% Loop over different substrate concentrations
for i = 1:length(C_s_range)
    C_s = C_s_range(i);
    
    % Define the Michaelis-Menten function
    v = @(P, S) (V_max * S) / (K_m + S);
    
    % Initial conditions
    P0 = 0; % Initial concentration of product (mg/mL)
    
    % Define the ODEs
    ode = @(t, Y) [R - v(Y(1), C_s); -v(Y(1), C_s)];
    
    % Solve the ODE using ode45
    [t, Y] = ode45(ode, tspan, [P0; C_s]);
    
    % Extract final concentration and reaction rate
    P_values(i) = Y(end, 1);
    ReactionRate_values(i) = V_max * C_s / (K_m + C_s);
end

% Specify a constant C_s value for the drug concentration plot
C_s_specified = 2.5; % ng/mL

% Calculate the drug concentration over time for the specified C_s
v_specified = @(P) (V_max * C_s_specified) / (K_m + C_s_specified);
ode_specified = @(t, P) R - v_specified(P);
[t_specified, P_specified] = ode45(ode_specified, tspan, P0);

% Plot the concentration vs. substrate concentration
figure;
subplot(3, 1, 1);
plot(C_s_range, P_values, 'b', 'LineWidth', 2);
xlabel('Substrate Concentration (ng/mL)');
ylabel('Isonicotinic Acid Hydrazide Concentration (mg/mL)');
title('Concentration of Isonicotinic Acid Hydrazide vs. Substrate Concentration');
grid on;

% Plot the reaction rate vs. substrate concentration
subplot(3, 1, 2);
plot(C_s_range, ReactionRate_values, 'r', 'LineWidth', 2);
xlabel('Substrate Concentration (ng/mL)');
ylabel('Reaction Rate (mg/mL/hr)');
title('Reaction Rate vs. Substrate Concentration');
grid on;

% Plot the concentration of the drug over time for the specified C_s
subplot(3, 1, 3);
plot(t_specified, P_specified, 'g', 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Isonicotinic Acid Hydrazide Concentration (mg/mL)');
title(['Concentration of Isonicotinic Acid Hydrazide Over Time for C_s = ', num2str(C_s_specified), ' ng/mL']);
grid on;