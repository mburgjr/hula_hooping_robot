% This is the main MATLAB script for Lab 5.
%
% You will need to modify the Mbed code and this script, but should not need to make any other changes.
%
%% SET YOUR INPUTS HERE

% Hoop
R_hoop = 0.24; % m

% Person
R_person = 0.03; % m

% Spiral input
% NOTE: Keep upper / lower equal for simple ellipse
a_upper = R_hoop/5; % m     Initial radius on x-axis
a_lower = R_hoop/5; % m     Final radius on x-axis
b_upper = R_hoop/5; % m     Initial radius on y-axis
b_lower = R_hoop/5; % m     Final radius on y-axis

a_spiral = [a_upper a_lower];
b_spiral = [b_upper b_lower];
v_spiral = 3; % m/s         Traversal speed
        
% Initial values of gen coords
h_init = 0;
phi_init = 0;

% Total experiment time is buffer,trajectory,buffer
pre_buffer_time   = 1;
traj_time         = 3;
post_buffer_time  = 1;

% Gains for impedance controller
% If a gain is not being used in your Mbed code, set it to zero
% For joint space control, use K_xx for K1, K_yy for K2, D_xx for D1, D_yy for D2
gains.K_xx = 200.0;
gains.K_yy = 200.0;
gains.K_xy = 10.0;

gains.D_xx = 15.0;
gains.D_yy = 15.0;
gains.D_xy = 2.5;

% Maximum duty cycle commanded by controller (should always be <=1.0)
duty_max   = 0.4;

%% Run Experiment
[output_data] = Experiment_trajectory( a_spiral, b_spiral, v_spiral,...
                                       h_init, phi_init,...
                                       traj_time, pre_buffer_time, post_buffer_time,...
                                       gains, duty_max);

%% Extract data
t = output_data(:,1);
x = -output_data(:,12); % actual contact position in X (negative due to direction motors are mounted)
y = output_data(:,13); % actual contact position in Y
   
xdes = -output_data(:,16); % desired contact position in X (negative due to direction motors are mounted)
ydes = output_data(:,17); % desired contact position in Y

%% Plot contact point vs desired
figure(3); clf;
subplot(211); hold on
plot(t,xdes,'r-'); plot(t,x);
xlabel('Time (s)'); ylabel('X (m)'); legend({'Desired','Actual'});

subplot(212); hold on
plot(t,ydes,'r-'); plot(t,y);
xlabel('Time (s)'); ylabel('Y (m)'); legend({'Desired','Actual'});

figure(4); clf; hold on
plot(xdes,ydes,'r-'); plot(x,y,'k');
xlabel('X (m)'); ylabel('Y (m)'); legend({'Desired','Actual'});
