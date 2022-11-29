function output_data = Experiment_trajectory( a_spiral, b_spiral, v_spiral, h_init, phi_init, traj_time, pre_buffer_time, post_buffer_time, gains, duty_max)
    
    % Figure for plotting motor data
    figure(1);  clf;       
    a1 = subplot(421);
    h1 = plot([0],[0]);
    h1.XData = []; h1.YData = [];
    ylabel('Phi (rad)');
    
    a2 = subplot(423);
    h2 = plot([0],[0]);
    h2.XData = []; h2.YData = [];
    ylabel('Phi dot (rad/s)');
    
    a3 = subplot(425);
    h3 = plot([0],[0]);
    h3.XData = []; h3.YData = [];
    ylabel('Current 1 (A)');
    hold on;
    subplot(425);
    h4 = plot([0],[0],'r');
    h4.XData = []; h4.YData = [];
    hold off;
    
    a4 = subplot(427);
    h5 = plot([0],[0]);
    h5.XData = []; h5.YData = [];
    ylabel('Duty Cycle 1');
    
    a5 = subplot(422);
    h21 = plot([0],[0]);
    h21.XData = []; h21.YData = [];
    ylabel('Theta (rad)');
    
    a6 = subplot(424);
    h22 = plot([0],[0]);
    h22.XData = []; h22.YData = [];
    ylabel('Theta Dot (rad/s)');
    
    a7 = subplot(426);
    h23 = plot([0],[0]);
    h23.XData = []; h23.YData = [];
    ylabel('Current 2 (A)');
    hold on;
    subplot(426);
    h24 = plot([0],[0],'r');
    h24.XData = []; h24.YData = [];
    hold off;
    
    a8 = subplot(428);
    h25 = plot([0],[0]);
    h25.XData = []; h25.YData = [];
    ylabel('Duty Cycle 2');
    
    % Figure for plotting state of the leg
    figure(2)
    clf
    hold on
    axis equal
    axis([-.35 .35 -.35 .1]);
   
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_ellip = plot([0],[0],'g','LineWidth',1.5);
    
    h_foot= plot([0],[0],'Color',[0.7,0.7,0.7]);
    h_des = plot([0],[0],'--','Color',[0.5,0.5,0.5]);
    h_des.XData=[];
    h_des.YData=[];
    h_foot.XData=[];
    h_foot.YData=[];
    
    % Define leg length parameters
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    Nmot = 18.75;
    Ir = 0.0035/Nmot^2;

    p   = [l_OA l_OB l_AC l_DE];
    p   = [p m1 m2 m3 m4 I1 I2 I3 I4 Ir Nmot l_O_m1 l_B_m2 l_A_m3 l_C_m4]';
    
    % This function will get called any time there is new data from
    % the Nucleo board. Data comes in blocks, rather than one at a time.
    function my_callback(new_data)
        % Parse new data
        t = new_data(:,1);          % time
        pos1 = new_data(:,2);       % position
        vel1 = new_data(:,3);       % velocity
        cur1 = new_data(:,4);       % current
        dcur1 = new_data(:,5);      % desired current
        duty1 = new_data(:,6);      % command
        
        pos2 = new_data(:,7);       % position
        vel2 = new_data(:,8);       % velocity
        cur2 = new_data(:,9);       % current
        dcur2 = new_data(:,10);     % desired current
        duty2 = new_data(:,11);     % command
        
        x = -new_data(:,12);         % actual foot position (negative due to direction motors are mounted)
        y = new_data(:,13);         % actual foot position
        xdes = -new_data(:,16);      % desired foot position (negative due to direction motors are mounted)
        ydes = new_data(:,17);      % desired foot position         
        
        N = length(pos1);
        
        % Update motor data plots
        h1.XData(end+1:end+N) = t;   
        h1.YData(end+1:end+N) = -pos1; % switch sign on all plotted values due to direction motors are mounted
        h2.XData(end+1:end+N) = t;   
        h2.YData(end+1:end+N) = -vel1;
        h3.XData(end+1:end+N) = t;   
        h3.YData(end+1:end+N) = -cur1;
        h4.XData(end+1:end+N) = t;   
        h4.YData(end+1:end+N) = -dcur1;
        h5.XData(end+1:end+N) = t;   
        h5.YData(end+1:end+N) = -duty1;
        
        h21.XData(end+1:end+N) = t;   
        h21.YData(end+1:end+N) = -pos2;
        h22.XData(end+1:end+N) = t;   
        h22.YData(end+1:end+N) = -vel2;
        h23.XData(end+1:end+N) = t;   
        h23.YData(end+1:end+N) = -cur2;
        h24.XData(end+1:end+N) = t;   
        h24.YData(end+1:end+N) = -dcur2;
        h25.XData(end+1:end+N) = t;   
        h25.YData(end+1:end+N) = -duty2;
        
        h_foot.XData(end+1:end+N) = x;
        h_foot.YData(end+1:end+N) = y;
        h_des.XData(end+1:end+N) = xdes;
        h_des.YData(end+1:end+N) = ydes;
        
    end
    
    frdm_ip  = '192.168.1.100';     % FRDM board ip
    frdm_port= 11223;               % FRDM board port  
    params.callback = @my_callback; % callback function
    % params.timeout  = 2;            % end of experiment timeout
    
    % Parameters for tuning
    start_period                = pre_buffer_time;    % In seconds 
    end_period                  = post_buffer_time;   % In seconds
    
    K_xx                     = gains.K_xx; % Stiffness
    K_yy                     = gains.K_yy; % Stiffness
    K_xy                     = gains.K_xy; % Stiffness

    D_xx                     = gains.D_xx; % Damping
    D_yy                     = gains.D_yy; % Damping
    D_xy                     = gains.D_xy; % Damping
    
    % Specify inputs
    input = [start_period traj_time end_period];
    input = [input a_spiral(1) a_spiral(2) b_spiral(1) b_spiral(2) v_spiral];
    input = [input h_init phi_init];
    input = [input K_xx K_yy K_xy D_xx D_yy D_xy];
    input = [input duty_max];
    
    params.timeout  = (start_period+traj_time+end_period);  
    
    output_size = 19;    % number of outputs expected
    output_data = RunExperiment(frdm_ip,frdm_port,input,output_size,params);
    linkaxes([a1 a2 a3 a4],'x')
    
end