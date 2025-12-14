function y = ZeroDynamics()
clear; close all; clc;

%% Define Parameters from Spong
m1 = 1;           % Mass of link 1
m2 = 1;           % Mass of link 2 (NOTE: Spong uses m2=1, not 2!)
l1 = 1;           % Length of link 1
l2 = 2;           % Length of link 2
l_cl1 = 0.5;      % COM location for link 1
l_cl2 = 1;        % COM location for link 2
I1 = 0.083;       % MOI for link 1
I2 = 0.33;        % MOI for link 2
g = 9.8;          % Gravity (Spong uses 9.8)

params = [m1, m2, l1, l2, l_cl1, l_cl2, I1, I2, g];

% Zero Dynamics aka equation 23 from Spong paper
% (m2*l_cl2^2 + m2*l1*l_cl2*cos(q2) + I2)*ddq2 - m2*l1*l_cl2*sin(q2)*dq2^2 
%                                              - m2*l_cl2*g*sin(q2) = 0


fig = figure('Position', [100 100 800 600]);
fig.Color = 'w';

hold on;
box on;
grid off;

% Plotting trajectories at differening energy levels 

% Defining energy levels at center at q2=0 and saddles at q2 = +- pi

energy_levels_center = [0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.5, 12.0, 15.0];


energy_levels_saddle = [0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.5, 12.0, 15.0];

tspan = [0 20];
line_color = [0 0.4470 0.7410]; % Blue 
first_traj = true; % For legend

%Energy countours around center
for E_target = energy_levels_center
    % First find initial angle with values q2_0 and dq2_0=0 giving
    % specified energy
  
    q2_0 = initial_energy(E_target, 0, params);
    
    if ~isnan(q2_0) && abs(q2_0) < pi
        [t, z] = ode45(@(t,z) zero_dynamics_ode(t, z, params), tspan, [q2_0; 0]);
        
        %Trajectory plot
        if first_traj
            plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
                'DisplayName', 'Zero Dynamics Trajectories');
            first_traj = false;
        else
            plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end
end

%Energy contours around saddles
for E_target = energy_levels_saddle
    % Find initial condition near saddle
    q2_0 = initial_energy(E_target, pi, params);
    
    if ~isnan(q2_0) && q2_0 > pi/2 && q2_0 < 3*pi/2
        [t, z] = ode45(@(t,z) zero_dynamics_ode(t, z, params), tspan, [q2_0; 0]);
        
        % Plot trajectory around pi
        plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
    
    % Trajectory around -pi.
    q2_0_neg = -q2_0;
    if ~isnan(q2_0_neg) && q2_0_neg < -pi/2 && q2_0_neg > -3*pi/2
        [t, z] = ode45(@(t,z) zero_dynamics_ode(t, z, params), tspan, [q2_0_neg; 0]);
        
        plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
end

% Trajectories initial velocity for connecting regions between the centers and saddles
dq2_initials = [0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.5, 12.0, 15.0];
q2_starts = [-pi/2, 0, pi/2];

for q2_0 = q2_starts
    for dq2_0 = dq2_initials

        %Positive velocities
        [t, z] = ode45(@(t,z) zero_dynamics_ode(t, z, params), [0 15], [q2_0; dq2_0]);
        plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        
        % Negative velocities
        [t, z] = ode45(@(t,z) zero_dynamics_ode(t, z, params), [0 15], [q2_0; -dq2_0]);
        plot(z(:,1), z(:,2), 'Color', line_color, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
end

% Plotting Center, Saddles, and Equillibria
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2, ...
    'DisplayName', 'Center (q_2=0)');
plot(pi, 0, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, ...
    'DisplayName', 'Saddle (q_2=\pm\pi)');
plot(-pi, 0, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, ...
    'HandleVisibility', 'off');
plot(2*pi, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2, ...
    'HandleVisibility', 'off');
plot(-2*pi, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2, ...
    'HandleVisibility', 'off');

%Plot format
xlabel('q2 (rad)', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
ylabel('dq2 (rad/s)', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
title('Phase portrait of zero dynamics for q2', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

xlim([-2*pi, 2*pi]);
ylim([-20, 20]);


set(gca, 'XTick', [-2*pi, -pi, 0, pi, 2*pi]);
set(gca, 'XTickLabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
set(gca, 'YTick', [-20, -10, 0, 10, 20]);
set(gca, 'FontSize', 11);
set(gca, 'LineWidth', 1.5);
set(gca, 'Color', 'w');  
set(gca, 'XColor', 'k');  
set(gca, 'YColor', 'k'); 

%Grid
grid on;
set(gca, 'GridColor', [0.85 0.85 0.85]);
set(gca, 'GridAlpha', 0.4);

%Legend
leg = legend('Location', 'northeast', 'FontSize', 10);
set(leg, 'TextColor', 'White');


hold off;

%% Functions
function dzdt = zero_dynamics_ode(t, z, params)
    % State vector: z = [q2; dq2]
    % Returns: dzdt = [dq2; ddq2]
    
    q2 = z(1);
    dq2 = z(2);
    
    m2 = params(2);
    l1 = params(3);
    l_cl2 = params(6);
    I2 = params(8);
    g = params(9);
    
    % Inertia coefficient
    J = m2*l_cl2^2 + m2*l1*l_cl2*cos(q2) + I2;
    
    % Zero dynamics
    coriolis_term = m2*l1*l_cl2*sin(q2)*dq2^2;
    gravity_term = m2*l_cl2*g*sin(q2);
    
    ddq2 = (coriolis_term + gravity_term) / J;
    
    dzdt = [dq2; ddq2];
end


function q2_0 = initial_energy(E_target, center, params)
    % Find initial position q2_0 (with dq2_0=0) that gives target energy
    % center: 0 for orbits around q2=0, pi for orbits around q2=pi
    
    m2 = 1;
    l_cl2 = 1;
    g = 9.81;
    
    % At dq2=0, E = V = -m2*g*l_cl2*cos(q2_0) (assumption from Spong eq 37)
    % cos(q2_0) = -E_target / (m2*g*l_cl2)
    
    cos_q2 = -E_target / (m2*g*l_cl2);
    
    if abs(cos_q2) > 1
        q2_0 = NaN;
        return;
    end
    
   
    if abs(center) < pi/4  
        q2_0 = acos(cos_q2);
    else  
        q2_0 = pi - acos(-cos_q2);
        if center < 0
            q2_0 = -q2_0;
        end
    end
end
y =0;
end
