%% Parameters

% Hoop
m_hoop = 1; % kg
R_hoop = 0.24; % m
I_hoop = m_hoop*R_hoop^2;

% Person
R_person = 0.01; % m
mu = 0.1; % [/]

% Elliptical input
a = R_hoop/2; % m       Length of x-axis limits
b = R_hoop/2; % m       Length of y-axis limits
dth = 2*pi; % rad/sec   Traversal speed

%% Forward simulation
% From person to hoop

t_lim = [0; 3];
dt = 0.01;
t = t_lim(1):dt:t_lim(2);
N = length(t);

% Generate ellipse
ang = t*dth;
p_person = [a*cos(ang);     b*sin(ang)];            % [x, y] X N
v_person = [-a*sin(ang);    b*cos(ang)]*dth;        % [dx, dy] X N
a_person = [-a*cos(ang);    -b*sin(ang)]*dth*dth;   % [ddx, ddy] X N





% State

% [th, r, dth, dr] of hoop

% Start here:
% Fix contact point
% Choose force within friction cone




% Init hoop states
p_hoop = zeros([3 N]);      % [x, y, phi] x N
v_hoop = zeros([3 N]);      % [dx, dy, dphi] x N

% Set initial state of hoop
p_hoop(1:2, 1) = p_person(:,1) + [R_person - R_hoop; 0];

% Iterate for timesteps
for i = 1:N



    % TODO: This loop is broken
    % It should use person trajectory to update hoop positions
    % While maintaining contact

    % My idea was that the interaction point supplies some force
    % To accelerate the hoop in the same direction as the person
    % And, a torque is applied about the center of the hoop according to
    % this acceleration vector of the person

    % 1. I don't think this is dynamically correct
    % 2. Torque hard :(



    % Find contact point, assuming contact
    % HP = p_person(:,i) - p_hoop(1:2,i);
    % HP_hat = HP / norm(HP);
    % contact_point = R_hoop*HP_hat + p_hoop(1:2,i);

    % Force calc
    % F_c = m_hoop * a_person(:,i);


    phi = p_hoop(3, i);
    F_c = [10*cos(phi); 10*sin(phi)];



    % Moment calc
    tau_c = mu*10*R_hoop;


    
    % Update hoop state
    v_hoop(1:2, i+1) = dt*F_c/m_hoop + v_hoop(1:2,i);
    v_hoop(3, i+1) = dt*tau_c/I_hoop + v_hoop(3,i);
    p_hoop(:, i+1) = p_hoop(:,i) + dt*v_hoop(:,i+1);

end

%% Backward simulation
% From hoop to person

% TODO? Lol

%% Display

% Prepare plot handles
hold on
h_hoop = plot([0],[0],'LineWidth',2);
h_person = plot([0],[0],'LineWidth',2);

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

% For visualization purposes
sim_speed = 0.05;

% Step through and update animation
for i = 1:N
    set(h_title,'String',  sprintf('t=%.2f',t(i)) ); % update title
    
    set(h_person,'XData', person_circle_x + p_person(1,i));
    set(h_person,'YData', person_circle_y + p_person(2,i));
    
    set(h_hoop,'XData', hoop_circle_x + p_hoop(1,i));
    set(h_hoop,'YData', hoop_circle_y + p_hoop(2,i));

    pause(dt/sim_speed)
end







