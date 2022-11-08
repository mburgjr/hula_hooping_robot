%% Parameters

% Hoop
m_hoop = 1; % kg
R_hoop = 0.24; % m
I_hoop = m_hoop*R_hoop^2;

% Person
R_person = 0.01; % m
mu = 0.1; % [/]

% Elliptical input
R_traj_person = R_hoop/2;
R_traj_hoop = R_traj_person - R_person + R_hoop;
dth = 2*pi; % rad/sec   Traversal speed

%% Simulation
% From person to hoop

t_lim = [0; 3];
dt = 0.01;
t = t_lim(1):dt:t_lim(2);
N = length(t);

% Generate ellipse
ang = t*dth;
p_person = [cos(ang);     sin(ang)]*R_traj_person;            % [x, y] X N
v_person = [-sin(ang);    cos(ang)]*R_traj_person*dth;        % [dx, dy] X N
a_person = [-cos(ang);   -sin(ang)]*R_traj_person*dth*dth;   % [ddx, ddy] X N

% Init hoop states
p_hoop = zeros([3 N]);      % [x, y, phi] x N
v_hoop = zeros([3 N]);      % [dx, dy, dphi] x N

% Set initial state of hoop
p_hoop(1:2, 1) = p_person(:,1) + [-R_person + R_hoop; 0];
p_hoop(3, 1) = pi;
v_hoop(2, 1) = p_hoop(1,1)*dth;
v_hoop(3, 1) = dth;

% Save contact point for plotting
p_contact = zeros([2 N]);

% Centripetal force required
T = 2*pi/dth;
F_centripetal = 4*pi^2*m_hoop*R_traj_hoop/(T^2);

% Iterate for timesteps
for i = 1:N
    
    % Apply centripetal force to inside of trajectory
    phi = p_hoop(3, i);
    F_c = F_centripetal*[cos(phi); sin(phi)];

    % Save contact point
    p_contact(:, i) = p_hoop(1:2,i) + [R_hoop*cos(phi); R_hoop*sin(phi)];

    % Moment calc
    tau_c = 0;

    % Update hoop state
    v_hoop(1:2, i+1) = dt*F_c/m_hoop + v_hoop(1:2,i);
    v_hoop(3, i+1) = dt*tau_c/I_hoop + v_hoop(3,i);
    p_hoop(:, i+1) = p_hoop(:,i) + dt*v_hoop(:,i);

end

%% Display

% Prepare plot handles
hold on
h_hoop = plot([0],[0],'b-','LineWidth',2);
h_person = plot([0],[0],'r-','LineWidth',2);
h_contact = plot([0],[0],'r.','MarkerSize',25);

xlabel('x'); ylabel('y');
h_title = title('t=0.0s');

axis equal
axis([-0.5 0.5 -0.5 0.5]);

% Precompute geometries
th = 0:pi/50:2*pi;
person_circle_x = R_person*cos(th);
person_circle_y = R_person*sin(th);
hoop_circle_x = R_hoop*cos(th);
hoop_circle_y = R_hoop*sin(th);

traj_circle_x = R_traj_person*cos(th);
traj_circle_y = R_traj_person*sin(th);
plot(traj_circle_x, traj_circle_y, 'k--')

traj_circle_x = R_traj_hoop*cos(th);
traj_circle_y = R_traj_hoop*sin(th);
plot(traj_circle_x, traj_circle_y, 'k--')

% For visualization purposes
sim_speed = 0.05;

% Step through and update animation
for i = 1:N
    set(h_title,'String',  sprintf('t=%.2f',t(i)) ); % update title
    
    % set(h_person,'XData', person_circle_x + p_person(1,i));
    % set(h_person,'YData', person_circle_y + p_person(2,i));
    
    set(h_hoop,'XData', hoop_circle_x + p_hoop(1,i));
    set(h_hoop,'YData', hoop_circle_y + p_hoop(2,i));

    set(h_contact,'XData', p_contact(1,i));
    set(h_contact,'YData', p_contact(2,i));

    pause(dt/sim_speed)
end







