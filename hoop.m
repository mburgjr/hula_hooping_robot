%% Parameters

% Hoop
m_hoop = 0.5; % kg
R_hoop = 0.24; % m
I_hoop = m_hoop*R_hoop^2;

% Person
R_person = 0.03; % m
mu = 0.1; % [/]

% Simulation
t_lim = [0; 3];
dt = 0.001;
t = t_lim(1):dt:t_lim(2);
N = length(t);

% Spiral input
% NOTE: Keep upper / lower equal for simple ellipse
a_upper = R_hoop; % m       Initial radius on x-axis
a_lower = R_hoop/10; % m    Final radius on x-axis
b_upper = R_hoop; % m       Initial radius on y-axis
b_lower = R_hoop/10; % m     Final radius on y-axis
v_spiral = 5; % m/s         Traversal speed

%% Generate spiral trajectory of person

a = [a_upper a_lower]; % Change of x-axis radius from start to end
b = [b_upper b_lower]; % Change of y-axis radius from start to end

% Calculate radius change over axes
dr_a = (a(2)-a(1)) / t_lim(2);
a_radius = a(1) + dr_a*t;
dr_b = (b(2)-b(1)) / t_lim(2);
b_radius = b(1) + dr_b*t;

% Calculate angle change
r_spiral = (a_radius.^2 + b_radius.^2).^0.5;
dth_spiral = v_spiral*(r_spiral.^-1);
th_spiral = zeros([1 N+1]);

p_person = zeros([2 N]);    % [x, y] x N
v_person = zeros([2 N]);    % [dx, dy] x N
a_person = zeros([2 N]);    % [ddx, ddy] x N

% Iteratively calculate spiral trajectory
for i = 1:N
    
    p_person(:,i) = [ a_radius(i)*cos(th_spiral(i)) ;...
                      b_radius(i)*sin(th_spiral(i)) ];
    v_person(:,i) = [ dr_a*cos(th_spiral(i)) - a_radius(i)*sin(th_spiral(i))*dth_spiral(i) ;...
                      dr_b*sin(th_spiral(i)) + b_radius(i)*cos(th_spiral(i))*dth_spiral(i) ];
    a_person(:,i) = [ dr_a*sin(th_spiral(i))*dth_spiral(i) - dr_a*sin(th_spiral(i))*dth_spiral(i) - a_radius(i)*cos(th_spiral(i))*dth_spiral(i)^2 ;...
                      dr_b*cos(th_spiral(i))*dth_spiral(i) + dr_b*cos(th_spiral(i))*dth_spiral(i) - b_radius(i)*sin(th_spiral(i))*dth_spiral(i)^2 ];

    % Update angle for next step
    th_spiral(:, i+1) = th_spiral(i) + dt*dth_spiral(i);
end

%% Forward simulation
% From person to hoop

% Init hoop states
p_hoop = zeros([3 N]);      % [x, y, phi] x N
v_hoop = zeros([3 N]);      % [dx, dy, dphi] x N

% Set initial state of hoop
p_hoop(1:2, 1) = p_person(:,1) + [-R_person + R_hoop; 0];
p_hoop(3, 1) = pi;

F_contact = zeros([2 N]);

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

        % Difference in momentum
        F_c = m_hoop*norm(v_plus - v_minus_hoop)*normal;

        % Calculate frictional torque
        tau_c = norm(F_c)*R_hoop*mu;

        % Tau is mu times normal component of F_c
        v_hoop(1:2, i+1) = F_c/m_hoop + v_hoop(1:2, i+1);
        v_hoop(3, i+1) = tau_c/I_hoop + v_hoop(3, i+1);

        % Save contact force for viewing
        F_contact(:,i) = F_c;
        
    end

    % Update hoop state
    p_hoop(:, i+1) = p_hoop(:,i) + dt*v_hoop(:,i+1);

end

%% Display

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
axis([-0.75 0.75 -0.75 0.75]);

% Precompute geometries
th = 0:pi/50:2*pi;
person_circle_x = R_person*cos(th);
person_circle_y = R_person*sin(th);
hoop_circle_x = R_hoop*cos(th);
hoop_circle_y = R_hoop*sin(th);

% For visualization purposes
sim_speed = 1;

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







