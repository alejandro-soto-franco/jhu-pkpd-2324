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
title('Concentration of Dextromethorphan over Time');
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
title('Logarithm of Concentration of Dextromethorphan over Time');
legend('Location', 'Best');