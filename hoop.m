%% Parameters

% Hoop
m_hoop = 0.5; % kg
R_hoop = 0.24; % m
I_hoop = m_hoop*R_hoop^2;

% Person
R_person = 0.03; % m
mu = 0.05; % [/]

% Simulation
t_lim = [0; 3];
dt = 0.001;
t = t_lim(1):dt:t_lim(2);
N = length(t);

% Spiral input (if not sweeping)
% NOTE: Keep upper / lower equal for simple ellipse
a_upper = R_hoop/2; % m       Initial radius on x-axis
a_lower = R_hoop/2; % m    Final radius on x-axis
b_upper = R_hoop/2; % m       Initial radius on y-axis
b_lower = R_hoop/2; % m     Final radius on y-axis
v_spiral = 3; % m/s         Traversal speed

% Run sweep
sweep = true;
R_sweep = R_hoop/10:R_hoop/20:R_hoop; % m
v_sweep = 0.75:0.125:5; % m/s

if sweep
    phase_diff_res = zeros([length(R_sweep) length(v_sweep)]);
    rise_time_res = zeros([length(R_sweep) length(v_sweep)]);
else
    phase_diff_res = zeros([1 1]);
    rise_time_res = zeros([1 1]);
end

%% Run simulation
for R_i = 1:size(phase_diff_res, 1)
    for v_i = 1:size(phase_diff_res, 2)
        
        if sweep
            % Unpack circular parameters
            a = [R_sweep(R_i) R_sweep(R_i)];
            b = [R_sweep(R_i) R_sweep(R_i)];
            v_spiral = v_sweep(v_i);
        else
            % Generate spiral trajectory of person
            a = [a_upper a_lower]; % Change of x-axis radius from start to end
            b = [b_upper b_lower]; % Change of y-axis radius from start to end
        end
        
        % Calculate radius change over axes
        dr_a = (a(2)-a(1)) / t_lim(2);
        radius_a = a(1) + dr_a*t;
        dr_b = (b(2)-b(1)) / t_lim(2);
        radius_b = b(1) + dr_b*t;
        
        % Calculate angle change
        r_spiral = (radius_a.^2 + radius_b.^2).^0.5;
        dth_spiral = v_spiral*(r_spiral.^-1);
        th_spiral = zeros([1 N+1]);
        
        p_person = zeros([2 N]);    % [x, y] x N
        v_person = zeros([2 N]);    % [dx, dy] x N
        a_person = zeros([2 N]);    % [ddx, ddy] x N
        ang_person = zeros([1 N]);  % [atan(y,x)] x N
        
        % Iteratively calculate spiral trajectory
        for i = 1:N
            r_a = radius_a(i);
            r_b = radius_b(i);
            th = th_spiral(i);
            dth = dth_spiral(i);
            
            p_person(:,i) = [ r_a*cos(th) ;...
                              r_b*sin(th) ];
            v_person(:,i) = [ dr_a*cos(th) - r_a*sin(th)*dth ;...
                              dr_b*sin(th) + r_b*cos(th)*dth ];
            a_person(:,i) = [ -dr_a*sin(th)*dth - dr_a*sin(th)*dth - r_a*cos(th)*dth^2 ;...
                              dr_b*cos(th)*dth + dr_b*cos(th)*dth - r_b*sin(th)*dth^2 ];
            ang_person(i) = atan2(p_person(2,i), p_person(1,i));
        
            % Update angle for next step
            th_spiral(:, i+1) = th + dt*dth;
        end
        
        % Contact dynamics
        K_c = 850; % N/m
        D_c = 85; % N*s/m
        
        % Init hoop states
        p_hoop = zeros([3 N]);      % [x, y, phi] x N
        v_hoop = zeros([3 N]);      % [dx, dy, dphi] x N
        ang_hoop = zeros([1 N]);    % [atan(y,x)] x N
        
        % Set initial state of hoop
        p_hoop(1:2, 1) = p_person(:,1) + [-R_person + R_hoop; 0];
        p_hoop(3, 1) = pi;
        ang_hoop(1) = 0;
        
        % Forward simulation (from person to hoop)
        F_contact = zeros([2 N]);
        ss_force_threshold = 150; % N
        phase_diff_log = zeros([1 N]);
        
        % Iterate for timesteps
        for i = 1:N-1
        
            % Carry over velocity
            v_hoop(:, i+1) = v_hoop(:,i);
        
            % Check if hoop and person are in contact
            dist_btwn = norm(p_person(:, i) - p_hoop(1:2, i));
            contact = (dist_btwn >= R_hoop - R_person) &&...
                        (dist_btwn <= R_hoop + R_person);
        
            % Apply force and torque to hoop if there's contact
            if contact
                % Planar collision
                v_plus = v_person(:, i+1);
                v_minus_hoop = v_hoop(1:2, i+1);
        
                % Collision force acts normal to the hoop
                normal = (p_person(:, i) - p_hoop(1:2, i)) / norm(p_person(:, i) - p_hoop(1:2, i));
        
                % Spring-damper collision model
                dx = dist_btwn - R_hoop + R_person;
                dv = norm(v_plus - v_minus_hoop);
                F_c = (K_c*dx + D_c*dv)*normal;
        
                % Tau is mu times normal component of F_c
                tau_c = norm(F_c)*R_hoop*mu;
                v_hoop(1:2, i+1) = dt*F_c/m_hoop + v_hoop(1:2, i+1);
                v_hoop(3, i+1) = dt*tau_c/I_hoop + v_hoop(3, i+1);
        
                % Save contact force for viewing
                F_contact(:,i) = F_c;
            end
        
            % Update hoop state
            p_hoop(:, i+1) = p_hoop(:,i) + dt*v_hoop(:,i+1);
            ang_hoop(i+1) = atan2(p_hoop(2,i+1), p_hoop(1,i+1));
        
            % Calculate and add to sum
            phase_diff_i = ang_hoop(i) - ang_person(i);
            phase_diff_i = mod((phase_diff_i + pi), 2*pi) - pi;
            phase_diff_log(i) = phase_diff_i;
        end

        % Determine when steady state is reached
        mavg_Fc = movmean(vecnorm(F_contact, 1),300);
        if sum(mavg_Fc > ss_force_threshold) > 0
            ss_index = find(mavg_Fc > ss_force_threshold,1) + 300;
            rise_time = dt*ss_index;
            phase_diff = mean(phase_diff_log(ss_index:end));
        else
            rise_time = inf;
            phase_diff = inf;
        end

        phase_diff_res(R_i, v_i) = phase_diff;
        rise_time_res(R_i, v_i) = rise_time;
    end
end

%% Display
if sweep == false
    
    watch_sim = true;
    if watch_sim
        % Plot trajectory
        hold on
        plot(p_person(1,:), p_person(2,:), '-', 'LineWidth', 1, 'Color', [0.65 0.65 0.65])
        
        % Prepare plot handles
        h_hoop = plot([0],[0],'b-','LineWidth',2);
        h_hoop_dot = plot([0],[0],'c.','MarkerSize',15);
        h_contact = plot([0],[0],'m-','LineWidth',2);
        h_person = fill([0],[0],'r-','LineStyle','none');
        
        xlabel('x'); ylabel('y');
        h_title = title('t=0.0s');
        
        axis equal
        axis([-0.65 0.65 -0.65 0.65]);
        
        % Precompute geometries
        th = 0:pi/50:2*pi;
        person_circle_x = R_person*cos(th);
        person_circle_y = R_person*sin(th);
        hoop_circle_x = R_hoop*cos(th);
        hoop_circle_y = R_hoop*sin(th);
        
        % For visualization purposes
        sim_speed = 10;
        
        % Step through and update animation
        for i = 1:N
            % Update time title
            set(h_title,'String',  sprintf('t=%.2f',t(i)) );
            
            % Plot bodies
            set(h_person,'XData', person_circle_x + p_person(1,i));
            set(h_person,'YData', person_circle_y + p_person(2,i));
            
            set(h_hoop,'XData', hoop_circle_x + p_hoop(1,i));
            set(h_hoop,'YData', hoop_circle_y + p_hoop(2,i));
        
            set(h_hoop_dot,'XData', [p_hoop(1,i) + R_hoop*cos(p_hoop(3,i))]);
            set(h_hoop_dot,'YData', [p_hoop(2,i) + R_hoop*sin(p_hoop(3,i))]);
        
            % Plot contact force
            normal = F_contact(:,i)/norm(F_contact(:,i));
            contact_vec_x = [0 0.15*normal(1)];
            contact_vec_y = [0 0.15*normal(2)];
            contact_point = p_hoop(1:2, i) + R_hoop*normal;
            set(h_contact,'XData', contact_vec_x + contact_point(1));
            set(h_contact,'YData', contact_vec_y + contact_point(2));
        
            % Pause sim between frames
            pause(dt/sim_speed)
        end
    
        hold off;
    end

    % Plot position
    figure;
    hold on;
    xlabel('Time (s)');
    ylabel('Position (m)');
    plot(dt*(1:N), p_person(1,:), 'r-');
    plot(dt*(1:N), p_person(2,:), 'm-');
    plot(dt*(1:N), p_hoop(1,:), 'b-');
    plot(dt*(1:N), p_hoop(2,:), 'c-');
    legend('Person X', 'Person Y', 'Hoop X', 'Hoop Y');
    hold off;

    % Plot velocity
    figure;
    hold on;
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    plot(dt*(1:N), v_person(1,:), 'r-');
    plot(dt*(1:N), v_person(2,:), 'm-');
    plot(dt*(1:N), v_hoop(1,:), 'b-');
    plot(dt*(1:N), v_hoop(2,:), 'c-');
    legend('Person X', 'Person Y', 'Hoop X', 'Hoop Y');
    hold off;

    % Plot angles
    figure;
    hold on;
    xlabel('Time (s)');
    ylabel('Ang (rad)');
    plot(dt*(1:N), ang_person(1,:), 'r-');
    plot(dt*(1:N), ang_hoop(1,:), 'b-');
    legend('Person', 'Hoop');
    hold off;

    % Plot contact force
    figure;
    hold on;
    xlabel('Time (s)');
    ylabel('Force (N)');
    plot(dt*(1:N), F_contact(1,:), 'r-');
    plot(dt*(1:N), F_contact(2,:), 'b-');
    legend('X', 'Y');
    hold off;

else
    % Plot sweep graphs
    figure, imagesc(R_sweep, v_sweep, rise_time_res');
    colormap([(0:0.01:1)', (1:-0.01:0)', 0.15*ones(101,1)]);
    colorbar;
    xlabel('Trajectory radius (m)');
    ylabel('Traversal speed (m/s)');
    title('Rise time of hoop');

    figure, imagesc(R_sweep, v_sweep, phase_diff_res');
    colormap([(0:0.01:1)', (1:-0.01:0)', 0.15*ones(101,1)]);
    colorbar;
    xlabel('Trajectory radius (m)');
    ylabel('Traversal speed (m/s)');
    title('Phase difference of hoop');
end






