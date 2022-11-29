#include "mbed.h"
#include "rtos.h"
#include "EthernetInterface.h"
#include "ExperimentServer.h"
#include "QEI.h"
#include "BezierCurve.h"
#include "MotorShield.h" 
#include "HardwareSetup.h"
#include "Matrix.h"
#include "MatrixMath.h"
#include "Servo.h"
#include "millis.h"

#define NUM_INPUTS 17
#define NUM_OUTPUTS 19
#define PI 3.14159265359

#define PULSE_TO_RAD (2.0f*3.14159f / 1200.0f)

// Initializations
Serial pc(USBTX, USBRX);    // USB Serial Terminal
ExperimentServer server;    // Object that lets us communicate with MATLAB
Timer t;                    // Timer to measure elapsed time of experiment
Servo sarrusservo(PB_3); //adjust PIN!

QEI encoderA(PE_9,PE_11, NC, 1200, QEI::X4_ENCODING);  // MOTOR A ENCODER (no index, 1200 counts/rev, Quadrature encoding) //Shaft motor encoder 



MotorShield motorShield(24000); //initialize the motor shield with a period of 12000 ticks or ~20kHZ
Ticker currentLoop;

Matrix MassMatrix(2,2);
Matrix Jacobian_(2,2);
Matrix JacobianT(2,2);
Matrix InverseMassMatrix(2,2);
Matrix temp_product(2,2);
Matrix Lambda(2,2);

// Variables for q1 (h)
float current1;
float current_des1 = 0;
float prev_current_des1 = 0;
float current_int1 = 0;
float h;
float velocity1;
float duty_cycle1;
float h_velocity;
float h_duty_cycle;
float h_init;
float angle1_init = 0; 

unsigned long timer;

// Variables for q2 (phi)
float current2;
float current_des2 = 0;
float prev_current_des2 = 0;
float current_int2 = 0;
float phi;
float phi_velocity;
float phi_duty_cycle;
float phi_init;

// Fixed kinematic parameters

// TODO: Change these to our values

// Lengths from point to point *
const float l_OA = 0.04064;
const float l_AB = 0.1016;
const float l_BC = l_AB;
const float l_CD = sqrt(pow(0.032,2) + pow(0.061,2)); 
const float l_DE = 0.096; 
const float l_EHoop=  0.022; // SUBJECT TO CHANGE// HOOP CONTACT POINT!
const float l_HoopG = 0.122 - l_EHoop; 
const float H = 0.4; // NEED TO MEASURE;

//jacobian calculation
// given x_F, y_F 
// phi = tan^-1(y_F/x_F) 
// now lets find h!
// find alpha (angle between hyp of F and l_FG)
//hyp = sqrt(x_F^2 +y_F^2);
// alpha = cos^-1(hyp/l_FG)
//find angle between leg1 and shaft
// beta = 90-alpha (degrees) 
// h1 = l_EG* cos(beta) 
//find l_shaftE 
//l_shaftE = sin(beta)*l_EG; 
// fin h2 
//h2 = sqrt(l_DE**2-l_shaftE**2)
// H = h1 +h2 




const float m1 =.0393 + .2;
const float m2 =.0368; 
const float m3 = .00783;
const float m4 = .0155;
const float I1 = 0.0000251;  //25.1 * 10^-6;
const float I2 = 0.0000535;  //53.5 * 10^-6;
const float I3 = 0.00000925; //9.25 * 10^-6;
const float I4 = 0.0000222;  //22.176 * 10^-6;
const float l_O_m1=0.032;
const float l_B_m2=0.0344; 
const float l_A_m3=0.0622;
const float l_C_m4=0.0610;
const float N = 18.75;
const float Ir = 0.0035/pow(N,2);

// Timing parameters
float current_control_period_us = 200.0f;     // 5kHz current control loop
float impedance_control_period_us = 1000.0f;  // 1kHz impedance control loop
float start_period, traj_time, end_period, a_upper, a_lower, b_upper, b_lower, v_spiral, traj_period;



// Control parameters
float current_Kp = 4.0f;         
float current_Ki = 0.4f;           
float current_int_max = 3.0f;       
float duty_max;      
float K_xx;
float K_yy;
float K_xy;
float D_xx;
float D_xy;
float D_yy;

// Model parameters
float supply_voltage = 12;     // motor supply voltage
float R = 2.0f;                // motor resistance
float k_t = 0.18f;             // motor torque constant
float nu = 0.0005;             // motor viscous friction

// Current control interrupt function
void CurrentLoop()
{
    // This loop sets the motor voltage commands using PI current controllers with feedforward terms.
    
    //use the motor shield as follows:
    //motorShield.motorAWrite(DUTY CYCLE, DIRECTION), DIRECTION = 0 is forward, DIRECTION =1 is backwards.
        
    current1 = -(((float(motorShield.readCurrentA())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c1 = current_des1 - current1;                                             // current errror
    current_int1 += err_c1;                                                             // integrate error
    current_int1 = fmaxf( fminf(current_int1, current_int_max), -current_int_max);      // anti-windup
    float ff1 = R*current_des1 + k_t*velocity1;                                         // feedforward terms
    duty_cycle1 = (ff1 + current_Kp*err_c1 + current_Ki*current_int1)/supply_voltage;   // PI current controller
    
    float absDuty1 = abs(duty_cycle1);
    if (absDuty1 > duty_max) {
        duty_cycle1 *= duty_max / absDuty1;
        absDuty1 = duty_max;
    }    
    if (duty_cycle1 < 0) { // backwards
        motorShield.motorAWrite(absDuty1, 1);
    } else { // forwards
        motorShield.motorAWrite(absDuty1, 0);
    }             
    prev_current_des1 = current_des1; 
    
    //current2     = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);       // measure current
//velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;                                  // measure velocity  
    //float err_c2 = current_des2 - current2;                                             // current error
    //current_int2 += err_c2;                                                             // integrate error
    //current_int2 = fmaxf( fminf(current_int2, current_int_max), -current_int_max);      // anti-windup   
    //float ff2 = R*current_des2 + k_t*velocity2;                                         // feedforward terms
    //duty_cycle2 = (ff2 + current_Kp*err_c2 + current_Ki*current_int2)/supply_voltage;   // PI current controller
    
    //float absDuty2 = abs(duty_cycle2);
//    if (absDuty2 > duty_max) {
//        duty_cycle2 *= duty_max / absDuty2;
//        absDuty2 = duty_max;
//    }    
//    if (duty_cycle2 < 0) { // backwards
//        motorShield.motorBWrite(absDuty2, 1);
//    } else { // forwards
//        motorShield.motorBWrite(absDuty2, 0);
//    }             
//    prev_current_des2 = current_des2; 
    
}

int main (void)
{
    
    // Object for 7th order Cartesian foot trajectory
    //BezierCurve rDesFoot_bez(2,BEZIER_ORDER_FOOT);
    
    // Link the terminal with our server and start it up
    server.attachTerminal(pc);
    server.init();
    
    // Continually get input from MATLAB and run experiments
    float input_params[NUM_INPUTS];
    pc.printf("%f",input_params[0]);
    
    while(1) {
        
        // If there are new inputs, this code will run
        if (server.getParams(input_params,NUM_INPUTS)) {
            
                        
            // Get inputs from MATLAB          
            start_period                = input_params[0];    // First buffer time, before trajectory (s)
            traj_time                   = input_params[1];    // Trajectory time/length (s)
            end_period                  = input_params[2];    // Second buffer time, after trajectory (s)
            
            a_upper                     = input_params[3];    // Initial radius on x-axis (m)
            a_lower                     = input_params[4];    // Final radius on x-axis (m)
            b_upper                     = input_params[5];    // Initial radius on y-axis (m)
            b_lower                     = input_params[6];    // Final radius on y-axis (m)
            v_spiral                    = input_params[7];    // Traversal speed (m/s)
    
            h_init                      = input_params[8];    // Initial height h (m)
            phi_init                    = input_params[9];    // Initial angle for shaft (rad)

            K_xx                        = input_params[10];    // Foot stiffness N/m
            K_yy                        = input_params[11];    // Foot stiffness N/m
            K_xy                        = input_params[12];    // Foot stiffness N/m
            D_xx                        = input_params[13];    // Foot damping N/(m/s)
            D_yy                        = input_params[14];    // Foot damping N/(m/s)
            D_xy                        = input_params[15];    // Foot damping N/(m/s)
            duty_max                    = input_params[16];    // Maximum duty factor

            // Unpack input spiral
            float dr_a = (a_lower - a_upper) / traj_time;
            float dr_b = (b_lower - b_upper) / traj_time;

            // Create trajectory to interpolate from
            int N = 500;
            float dt = traj_time / N;
            float spiral_t[N];
            float spiral_x[N];
            float spiral_y[N];
            float spiral_xdot[N];
            float spiral_ydot[N];

            float th = 0;
            for(int i = 0; i<N; i++) {
                float t = i*dt;
                float radius_a = a_upper + dr_a*t;
                float radius_b = b_upper + dr_a*t;
                float r_spiral = sqrt(pow(radius_a,2) + pow(radius_b,2));

                float dth_spiral = v_spiral* pow(r_spiral,-1);
                th += dth_spiral*dt;
            
                spiral_t[i] = t;
                spiral_x[i] = radius_a*cos(th);
                spiral_y[i] = radius_b*sin(th);
                spiral_xdot[i] = dr_a*cos(th) - radius_a*sin(th)*dth_spiral;
                spiral_ydot[1] = dr_b*sin(th) + radius_b*cos(th)*dth_spiral;
            };




            // NOTE: I changed angle1 => h and angle2 => phi
            
            



            // Attach current loop interrupt
            currentLoop.attach_us(CurrentLoop,current_control_period_us);
                        
            // Setup experiment
            t.reset();
            t.start();
            encoderA.reset();
        
            motorShield.motorAWrite(0, 0); //turn motor A off
            // motorShield.motorBWrite(0, 0); //turn motor B off
            sarrusservo.write(90); //sets servo to midpoint
            float angle2_prev = 90;
            timer = millis(); 

            // Run experiment
            while( t.read() < start_period + traj_period + end_period) { 
                 
                // Read encoders to get motor states
                float angle1 = encoderA.getPulses() *PULSE_TO_RAD + angle1_init;       
                float velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;

                //servo current angle 
                float angle2 = sarrusservo.read(); 
                dt = millis()-timer;
                timer = millis(); 
                //servo velocity
                float velocity2 = (angle2-angle2_prev)/dt; //Servo angular velocity but unsure how to get this 
                angle2_prev = angle2; 

                 
                // angle2 = encoderB.getPulses() * PULSE_TO_RAD + angle2_init;       
                // velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;           
                
                const float th1 = angle1; //phi
                const float th2 = angle2; //theta
                const float dth1= velocity1;
                const float dth2= velocity2;
 
                // Calculate the Jacobian 1 (from phi, theta to phi, h)

                float h = H - 2*l_AB*cos(th2);
                float dh = 2*l_AB*sin(th2)*dth2;
                float phi = th1;
                float dphi = dth1;

                float Jh_th1 = 0;
                float Jh_th2 = 2*l_AB*sin(th2);
                float Jphi_th1 = 1;
                float Jphi_th2 = 0;

                // Calculate the Jacobian 2 (from phi,h to xHoop, yHoop) 

                //These are really long and crazy - check the jacobian derivation in MATLAB code

                float Jx_h = (l_HoopG*sin(acos(pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2)/(2*h*l_HoopG)) + PI/2)*cos(phi)*(1/l_HoopG - (pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2)/(2*pow(h,2)*l_HoopG))))/(1 - (pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2)/pow((4*pow(h,2)*pow(l_HoopG,2)),0.5)));
                float Jx_phi= -l_HoopG*cos(acos((pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*h*l_HoopG)) + PI/2)*sin(phi);
                float Jy_h = (l_HoopG*sin(acos((pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*h*l_HoopG)) + PI/2)*sin(phi)*(1/l_HoopG - (pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*pow(h,2)*l_HoopG)))/pow(1 - pow(pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2), 2)/(4*pow(h,2)*pow(l_HoopG,2)), 0.5);
                float Jy_phi = l_HoopG*cos(acos((pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*h*l_HoopG)) + PI/2)*cos(phi);
                                
                //Calculate the total Jacobian (J2*J1)

                float Jx_th1 = Jh_th1*Jx_h + Jphi_th1*Jx_phi;
                float Jx_th2 = Jh_th2*Jx_h + Jphi_th2*Jx_phi;
                float Jy_th1 = Jh_th1*Jy_h + Jphi_th1*Jy_phi;
                float Jy_th2 = Jh_th2*Jy_h + Jphi_th2*Jy_phi;
                // Calculate the forward kinematics (position and velocity) // calculate xF and yF
                float xHoop = l_HoopG*cos(acos((pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*h*l_HoopG)) + PI/2)*cos(th1);
                float yHoop = l_HoopG*cos(acos((pow(h,2) - pow(l_DE,2) + pow(l_HoopG,2))/(2*h*l_HoopG)) + PI/2)*sin(th1);

                //These are really long and crazy - check the jacobian derivation in MATLAB code

                float dxHoop = 0;  //(dh*l_HoopG*sin(acos((h^2 - l_DE^2 + l_HoopG^2)/(2*h*l_HoopG)) + PI/2)*cos(th1)*(1/l_HoopG - (h^2 - l_DE^2 + l_HoopG^2)/(2*h^2*l_HoopG)))/(1 - (h^2 - l_DE^2 + l_HoopG^2)^2/(4*h^2*l_HoopG^2))^(1/2) - dth1*l_HoopG*cos(acos((h^2 - l_DE^2 + l_HoopG^2)/(2*h*l_HoopG)) + PI/2)*sin(phi);
                float dyHoop = 0; //dth1*l_HoopG*cos(acos((h^2 - l_DE^2 + l_HoopG^2)/(2*h*l_HoopG)) + PI/2)*cos(th1) + (dh*l_HoopG*sin(acos((h^2 - l_DE^2 + l_HoopG^2)/(2*h*l_HoopG)) + PI/2)*sin(th1)*(1/l_HoopG - (h^2 - l_DE^2 + l_HoopG^2)/(2*h^2*l_HoopG)))/(1 - (h^2 - l_DE^2 + l_HoopG^2)^2/(4*h^2*l_HoopG^2))^(1/2);
 
   

                // TO ADJUST 
                // Set gains based on buffer and traj times, then calculate desired x,y from Bezier trajectory at current time if necessary
                float teff  = 0;
                float vMult = 0;
                if( t < start_period) {
                    if (K_xx > 0 || K_yy > 0) {
                        K_xx = 100; 
                        K_yy = 100; 
                        D_xx = 5;  
                        D_yy = 5;  
                        K_xy = 0;
                        D_xy = 0;
                    }
                    teff = 0;
                }
                else if (t < start_period + traj_period)
                {
                    K_xx = input_params[5];  // Foot stiffness N/m
                    K_yy = input_params[6];  // Foot stiffness N/m
                    K_xy = input_params[7];  // Foot stiffness N/m
                    D_xx = input_params[8];  // Foot damping N/(m/s)
                    D_yy = input_params[9];  // Foot damping N/(m/s)
                    D_xy = input_params[10]; // Foot damping N/(m/s)
                    teff = (t-start_period);
                    vMult = 1;
                }
                else
                {
                    teff = traj_period;
                    vMult = 0;
                }
                
                // Get desired workspace point from spiral
                float rDesContact[2] , vDesContact[2];
                for (int i = 0; i<(N-1); i++) {
                    if (teff == spiral_t[i] | i == N-1) {
                        // Set exact value
                        rDesContact[0] = spiral_x[i];
                        rDesContact[1] = spiral_y[i];
                        vDesContact[0] = spiral_xdot[i];
                        vDesContact[1] = spiral_ydot[i];
                        break;
                    }
                    else if (teff >= spiral_t[i] & teff < spiral_t[i+1]) {
                        // Interpolate
                        float delta = (teff - spiral_t[i])/(spiral_t[i+1] - spiral_t[i]);
                        rDesContact[0] = spiral_x[i] + delta*(spiral_x[i+1] - spiral_x[i]);
                        rDesContact[1] = spiral_y[i] + delta*(spiral_y[i+1] - spiral_y[i]);
                        vDesContact[0] = spiral_xdot[i] + delta*(spiral_xdot[i+1] - spiral_xdot[i]);
                        vDesContact[1] = spiral_ydot[i] + delta*(spiral_ydot[i+1] - spiral_ydot[i]);
                        break;
                    };

                };
                



                // ADJUSTED
                // Calculate the inverse kinematics (joint positions and velocities) for desired joint angles              
                float xHoop_inv = -rDesContact[0];
                float yHoop_inv = rDesContact[1];                
                float r = sqrt(pow(xHoop_inv,2) + pow(yHoop_inv,2) );
                float gamma = abs(acos(r/l_HoopG)); 
                float alpha = 3.14159f-gamma; 
                float h_des1 = cos(alpha)*(l_EHoop + l_HoopG); 
                float straight_edge = sin(alpha)*(l_EHoop + l_HoopG); 
                float h_des2 = sqrt(pow(straight_edge,2) + pow(l_DE,2)); 
                float h_des = h_des1 + h_des2; 
                float th2_des = acos((H-h_des)/(2*l_AB)); 
                float th1_des = -(3.14159f/2.0f) + atan2(yHoop_inv,xHoop_inv); 
                
                float dd = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                //float dth1_des = (1.0f/dd) * (  Jy_th2*vDesHoop[0] - Jx_th2*vDesHoop[1] );
                //float dth2_des = (1.0f/dd) * ( -Jy_th1*vDesHoop[0] + Jx_th1*vDesHoop[1] );
        
                // Calculate error variables
                float e_x = rDesContact[0] - xHoop;
                float e_y = rDesContact[1] - yHoop;
                float de_x = vDesContact[0] - dxHoop;
                float de_y = vDesContact[1] - dyHoop;
        
                // Calculate virtual force on foot
                float fx = K_xx*e_x + K_xy*e_y + D_xx*de_x + D_xy*de_y;
                float fy = K_xy*e_x + K_yy*e_y + D_xy*de_x + D_yy*de_y;
                
                current_des1 = (Jx_th1*fx + Jy_th1*fy)/k_t;
                current_des2 = (Jy_th2*fy + Jx_th2*fx)/k_t; 

                // Form output to send to MATLAB     
                float output_data[NUM_OUTPUTS];
                // current time
                output_data[0] = t.read();
                // motor 1 state
                output_data[1] = h;
                output_data[2] = h_velocity;  
                output_data[3] = current1;
                output_data[4] = current_des1;
                output_data[5] = h_duty_cycle;
                // motor 2 state
                output_data[6] = phi;
                output_data[7] = phi_velocity;
                output_data[8] = current2;
                output_data[9] = current_des2;
                output_data[10]= phi_duty_cycle;
                // foot state
                output_data[11] = xHoop;
                output_data[12] = yHoop;
                output_data[13] = dxHoop;
                output_data[14] = dyHoop;
                output_data[15] = rDesContact[0];
                output_data[16] = rDesContact[1];
                output_data[17] = vDesContact[0];
                output_data[18] = vDesContact[1];
                
                // Send data to MATLAB
                server.sendData(output_data,NUM_OUTPUTS);

                wait_us(impedance_control_period_us);   
            }
            
            // Cleanup after experiment
            server.setExperimentComplete();
            currentLoop.detach();
            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
        
        } // end if
        
    } // end while
    
} // end main