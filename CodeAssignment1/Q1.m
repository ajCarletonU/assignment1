%ELEC4700 Winter 2022
%ELEC4700 Winter 2022
%Assignment 1

%Alina Jacobson     101055071

set(0,'DefaultFigureWindowStyle','docked')
close all;

%-------------------------------------------------------------------------
%Q1
%-------------------------------------------------------------------------

%define constant variables 
k_B = 1.38064852e-23;               % universal constant (Boltzman )
mass_o = 9.10938356e-31;            % in kg is the rest mass
m_n = 0.26*mass_o;                  % in kg is effective mass of electrons
tempK = 300;                         % in Kelvin (temperature at room)



%program that will model the random motion of the electrons. 
particle_count = 30;      %number of particles
particle_toPlot = 10;      %number of particles to plot in simulation
numItr =20;
Xlength = 2e-7;            %horizontal - based on figure2 in assignment1 given     
Yheight = 1e-7;            %vertical   - based on figure2 in assignment1 given  


%The thermal velocity vth - Assume T = 300 K. 
V_thermal = sqrt(2*k_B*tempK/m_n);       % in m/s 

%spacial step be smaller than 1/100 of the region which is 2e-7, 
% 1/100 *2e-7 = 2e-9 / Vth
dt_step_size = Yheight/(V_thermal*100);    %smallest step allowed
time_step = 1000;                       %Simulate for nominally 1000 timesteps.                        
max_time = dt_step_size*time_step;


%The mean free Path 
t_mn = 0.2e-12;                     %mean time between collision
mean_fp = V_thermal*t_mn;           %mean free path = 3.7404e-8


%Assign each particle a random location in the x−y plane within the region 
%defined by the extent of the Silicom
xcord(1:particle_count) = rand(particle_count,1)*Xlength;             % on the length of the region
ycord(1:particle_count) = rand(particle_count,1)*Yheight;             % on the height of the region


%Assign each initial velocity given by vth 
%but give each one a random direction.
%set initial velocities of the electrons at random
%find the total velocity
x_vel(1:particle_count) = V_thermal*cos(2*pi*randn(particle_count,1)); 
y_vel(1:particle_count) = V_thermal*sin(2*pi*randn(particle_count,1)); 


x_cordprev = zeros(particle_count);                     %use the previos instances of the trajectory of the electon position
y_cordprev = zeros(particle_count);  

%before generating values in the loop set all values to back to zero
 
tempMaterial = 0;                     
Vel_bounds_Avg = 0;                    
total_temp = 0;
tempMaterial_prev = 0;
count=0;

while count < max_time 
    
    x_cordprev = xcord;
    y_cordprev = ycord;
    
    xcord(1:particle_count) = xcord(1:particle_count) + (dt_step_size .* x_vel);  
    ycord(1:particle_count) = ycord(1:particle_count) + (dt_step_size .* y_vel);
    
        
    %loop through to get calculed value output up to a max threshold
    %loop through to check x and y boundries in the region
    for i=1:particle_count
      
        %x boundary
       if xcord(i) >= Xlength
           x_cordprev(i) = x_cordprev(i)-Xlength;
           xcord(i) = xcord(i)-Xlength;       
       end
       if xcord(i) <= 0
           x_cordprev(i) = x_cordprev(i) + Xlength;
           xcord(i) =xcord(i) +Xlength ;      
       end
       
       %y boundary
       if ycord(i) >=Yheight
           y_vel(i) = -y_vel(i);
       end
       if ycord(i)<= 0
           y_vel(i) = -y_vel(i);
       end
   
 
    end
    
          %Store valvue of the previos and current positions- for plotting
%         x = [x_cordprev(i) xcord(i)];   
%         y = [y_cordprev(i) ycord(i)]; 

       %do the temperature calc  an check bounds
       Vel_bounds_Avg = mean(sqrt(x_vel.^2 + y_vel.^2));       
       total_temp = (m_n*Vel_bounds_Avg.^2)/(2*k_B);
     
        figure(1) 
        hold on
        subplot(2,1,1)
        axis([0 Xlength 0 Yheight]);
        title(['Q1 Electron Modelling   (' ,num2str(particle_toPlot), ' electrons)']);
        xlabel('xVelocity (nm)');
        ylabel('yVeloncity (nm)');
        
   
         
        for z=1:particle_toPlot
           %display in subplot 1 the moving electrons             
           colors = hsv(particle_toPlot);          %set to have randoms colors of the stream
           plot([x_cordprev(z) xcord(z)],[y_cordprev(z) ycord(z)], '-', 'markers', 1, 'color', colors(z,:), 'MarkerFaceColor', colors(z,:));
           
        end
        
        
       Temp = [tempMaterial_prev total_temp];
       Time = [(count-dt_step_size) count];
       hold on  
       subplot(2,1,2)
 
        plot(Time,Temp,'r-')
        %axis([0 1000 0 500]);
        title(['Average Temperature:', num2str(total_temp), 'K']);
        xlabel('xTime(s)');
        ylabel('yTemperature (K)'); 
        
        
        tempMaterial_prev =total_temp; %set the prev temp as the new temp
        total_temp =0;                  %then clear for new value
     
       pause(0.0001); 
       count=count+dt_step_size;        %increment the count one step     
        
     

end   
 

