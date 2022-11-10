%% Parameters

% Hoop
m_hoop = 0.25; % kg
R_hoop = 0.24; % m
I_hoop = m_hoop*R_hoop^2;

% Person
R_person = 0.03; % m
mu = 0.1; % [/]

% Elliptical input
a = R_hoop/2; % m       Length of x-axis limits
b = R_hoop/2; % m       Length of y-axis limits
dth = 4*pi; % rad/sec   Traversal speed

%% Forward simulation
% From person to hoop

t_lim = [0; 3];
dt = 0.001;
t = t_lim(1):dt:t_lim(2);
N = length(t);

% Generate ellipse
ang = t*dth;
p_person = [a*cos(ang);     b*sin(ang)];            % [x, y] X N
v_person = [-a*sin(ang);    b*cos(ang)]*dth;        % [dx, dy] X N
a_person = [-a*cos(ang);    -b*sin(ang)]*dth*dth;   % [ddx, ddy] X N

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

        % TODO: Implement torque
        tau_c = 0;

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
sim_speed = 0.1;

% Step through and update animation
for i = 1:N
    set(h_title,'String',  sprintf('t=%.2f',t(i)) ); % update title
    
    set(h_person,'XData', person_circle_x + p_person(1,i));
    set(h_person,'YData', person_circle_y + p_person(2,i));
    
    set(h_hoop,'XData', hoop_circle_x + p_hoop(1,i));
    set(h_hoop,'YData', hoop_circle_y + p_hoop(2,i));

    contact_hat = 0.1*F_contact(:,i)/norm(F_contact(:,i));
    contact_vec_x = [0 contact_hat(1)];
    contact_vec_y = [0 contact_hat(2)];
    normal = (p_person(:, i) - p_hoop(1:2, i)) / norm(p_person(:, i) - p_hoop(1:2, i));
    contact_point = p_hoop(1:2, i) + R_hoop*normal;
    set(h_contact,'XData', contact_vec_x + contact_point(1));
    set(h_contact,'YData', contact_vec_y + contact_point(2));

    pause(dt/sim_speed)
end







