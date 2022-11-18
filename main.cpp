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

#define NUM_INPUTS 17
#define NUM_OUTPUTS 19

#define PULSE_TO_RAD (2.0f*3.14159f / 1200.0f)

// Initializations
Serial pc(USBTX, USBRX);    // USB Serial Terminal
ExperimentServer server;    // Object that lets us communicate with MATLAB
Timer t;                    // Timer to measure elapsed time of experiment

QEI encoderA(PE_9,PE_11, NC, 1200, QEI::X4_ENCODING);  // MOTOR A ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderB(PA_5, PB_3, NC, 1200, QEI::X4_ENCODING);  // MOTOR B ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderC(PC_6, PC_7, NC, 1200, QEI::X4_ENCODING);  // MOTOR C ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderD(PD_12, PD_13, NC, 1200, QEI::X4_ENCODING);// MOTOR D ENCODER (no index, 1200 counts/rev, Quadrature encoding)

MotorShield motorShield(24000); //initialize the motor shield with a period of 12000 ticks or ~20kHZ
Ticker currentLoop;

Matrix MassMatrix(2,2);
Matrix Jacobian(2,2);
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
float h_velocity;
float h_duty_cycle;
float h_init;

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

const float l_OA=.011; 
const float l_OB=.042; 
const float l_AC=.096; 
const float l_DE=.091;
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
float start_period, traj_time, end_period;

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
    
    current2     = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);       // measure current
    velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;                                  // measure velocity  
    float err_c2 = current_des2 - current2;                                             // current error
    current_int2 += err_c2;                                                             // integrate error
    current_int2 = fmaxf( fminf(current_int2, current_int_max), -current_int_max);      // anti-windup   
    float ff2 = R*current_des2 + k_t*velocity2;                                         // feedforward terms
    duty_cycle2 = (ff2 + current_Kp*err_c2 + current_Ki*current_int2)/supply_voltage;   // PI current controller
    
    float absDuty2 = abs(duty_cycle2);
    if (absDuty2 > duty_max) {
        duty_cycle2 *= duty_max / absDuty2;
        absDuty2 = duty_max;
    }    
    if (duty_cycle2 < 0) { // backwards
        motorShield.motorBWrite(absDuty2, 1);
    } else { // forwards
        motorShield.motorBWrite(absDuty2, 0);
    }             
    prev_current_des2 = current_des2; 
    
}

int main (void)
{
    
    // Object for 7th order Cartesian foot trajectory
    BezierCurve rDesFoot_bez(2,BEZIER_ORDER_FOOT);
    
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
            float N = 500;
            float dt = traj_time / N;
            float spiral_t[N], spiral_x[N], spiral_y[N], spiral_xdot[N], spiral_ydot[N];

            float th = 0;
            for(int i = 0; i<N; i++) {
                float t = i*dt;
                float radius_a = a_upper + dr_a*t;
                float radius_b = b_upper + dr_a*t;
                float r_spiral = (radius_a^2 + radius_b^2)^0.5;

                float dth_spiral = v_spiral*(r_spiral^-1);
                th += dth_spiral*dt;
            
                spiral_t[i] = t;
                spiral_x[i] = radius_a*cos(th);
                spiral_y[i] = radius_b*sin(th);
                spiral_xdot[i] = dr_a*cos(th) - radius_a*sin(th)*dth;
                spiral_ydot[1] = dr_b*sin(th) + radius_b*cos(th)*dth;
            };








            // NOTE: I changed angle1 => h and angle2 => phi
            
            



            // Attach current loop interrupt
            currentLoop.attach_us(CurrentLoop,current_control_period_us);
                        
            // Setup experiment
            t.reset();
            t.start();
            encoderA.reset();
            encoderB.reset();
            encoderC.reset();
            encoderD.reset();

            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
                         
            // Run experiment
            while( t.read() < start_period + traj_period + end_period) { 
                 
                // Read encoders to get motor states
                angle1 = encoderA.getPulses() *PULSE_TO_RAD + angle1_init;       
                velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;
                 
                angle2 = encoderB.getPulses() * PULSE_TO_RAD + angle2_init;       
                velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;           
                
                const float th1 = angle1;
                const float th2 = angle2;
                const float dth1= velocity1;
                const float dth2= velocity2;
 
                // Calculate the Jacobian
                float Jx_th1 = l_AC*cos(th1 + th2) + l_DE*cos(th1) + l_OB*cos(th1);
                float Jx_th2 = l_AC*cos(th1 + th2);
                float Jy_th1 = l_AC*sin(th1 + th2) + l_DE*sin(th1) + l_OB*sin(th1);
                float Jy_th2 = l_AC*sin(th1 + th2);
                                
                // Calculate the forward kinematics (position and velocity)
                float xFoot = l_AC*sin(th1 + th2) + l_DE*sin(th1) + l_OB*sin(th1);
                float yFoot = -l_AC*cos(th1 + th2) - l_DE*cos(th1) - l_OB*cos(th1);
                float dxFoot = dth1*(l_AC*cos(th1 + th2) + l_DE*cos(th1) + l_OB*cos(th1)) + dth2*l_AC*cos(th1 + th2);
                float dyFoot = dth1*(l_AC*sin(th1 + th2) + l_DE*sin(th1) + l_OB*sin(th1)) + dth2*l_AC*sin(th1 + th2);      

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
                




                // TODO: Interpolate values from unpacked spiral above to get desired points
                float rDesContact[2] , vDesContact[2];
                rDesContact[0] = ;
                rDesContact[1] = ;
                vDesContact[0] = ;
                vDesContact[1] = ;
                




                // Calculate the inverse kinematics (joint positions and velocities) for desired joint angles              
                float xFoot_inv = -rDesFoot[0];
                float yFoot_inv = rDesFoot[1];                
                float l_OE = sqrt( (pow(xFoot_inv,2) + pow(yFoot_inv,2)) );
                float alpha = abs(acos( (pow(l_OE,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th2_des = -(3.14159f - alpha); 
                float th1_des = -((3.14159f/2.0f) + atan2(yFoot_inv,xFoot_inv) - abs(asin( (l_AC/l_OE)*sin(alpha) )));
                
                float dd = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                float dth1_des = (1.0f/dd) * (  Jy_th2*vDesFoot[0] - Jx_th2*vDesFoot[1] );
                float dth2_des = (1.0f/dd) * ( -Jy_th1*vDesFoot[0] + Jx_th1*vDesFoot[1] );
        
                // Calculate error variables
                float e_x = rDesFoot[0] - xFoot;
                float e_y = rDesFoot[1] - yFoot;
                float de_x = vDesFoot[0] - dxFoot;
                float de_y = vDesFoot[1] - dyFoot;
        
                // Calculate virtual force on foot
                float fx = K_xx*e_x + K_xy*e_y + D_xx*de_x + D_xy*de_y;
                float fy = K_xy*e_x + K_yy*e_y + D_xy*de_x + D_yy*de_y;
                
                // Calculate mass matrix elements
                float M11 = I1 + I2 + I3 + I4 + Ir + Ir*pow(N,2) + pow(l_AC,2)*m4 + pow(l_A_m3,2)*m3 + pow(l_B_m2,2)*m2 + pow(l_C_m4,2)*m4 + pow(l_OA,2)*m3 + pow(l_OB,2)*m2 + pow(l_OA,2)*m4 + pow(l_O_m1,2)*m1 + 2*l_C_m4*l_OA*m4 + 2*l_AC*l_C_m4*m4*cos(th2) + 2*l_AC*l_OA*m4*cos(th2) + 2*l_A_m3*l_OA*m3*cos(th2) + 2*l_B_m2*l_OB*m2*cos(th2); 
                float M12 = I2 + I3 + pow(l_AC,2)*m4 + pow(l_A_m3,2)*m3 + pow(l_B_m2,2)*m2 + Ir*N + l_AC*l_C_m4*m4*cos(th2) + l_AC*l_OA*m4*cos(th2) + l_A_m3*l_OA*m3*cos(th2) + l_B_m2*l_OB*m2*cos(th2); 
                float M22 = Ir*pow(N,2) + m4*pow(l_AC,2) + m3*pow(l_A_m3,2) + m2*pow(l_B_m2,2) + I2 + I3;
                
                
                // Populate mass matrix
                MassMatrix.Clear();
                MassMatrix << M11 << M12
                           << M12 << M22;
                
                // Populate Jacobian matrix
                Jacobian.Clear();
                Jacobian << Jx_th1 << Jx_th2
                         << Jy_th1 << Jy_th2;
                
                // Once you have copied the elements of the mass matrix, uncomment the following section
                
                // Calculate Lambda matrix
                JacobianT = MatrixMath::Transpose(Jacobian);
                InverseMassMatrix = MatrixMath::Inv(MassMatrix);
                temp_product = Jacobian*InverseMassMatrix*JacobianT;
                Lambda = MatrixMath::Inv(temp_product); 
                
                // Pull elements of Lambda matrix
                float L11 = Lambda.getNumber(1,1);
                float L12 = Lambda.getNumber(1,2);
                float L21 = Lambda.getNumber(2,1);
                float L22 = Lambda.getNumber(2,2);               
                                
                                
                                
                // Set desired currents
                current_des1 = (fx*(Jx_th1*L11 + Jy_th1*L21) + fy*(Jx_th1*L12 + Jy_th1*L22))/k_t;          
                current_des2 = (fx*(Jx_th2*L11 + Jy_th2*L21) + fy*(Jx_th2*L12 + Jy_th2*L22))/k_t;   







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
                output_data[11] = xContact;
                output_data[12] = yContact;
                output_data[13] = dxContact;
                output_data[14] = dyContact;
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
