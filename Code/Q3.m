%ELEC4700 Winter 2022
%Assignment 1

%Alina Jacobson     101055071

set(0,'DefaultFigureWindowStyle','docked')
close all;

%-------------------------------------------------------------------------
%Q3
%-------------------------------------------------------------------------
%define constant variables 
k_B = 1.38064852e-23;               % universal constant (Boltzman )
mass_o = 9.10938356e-31;            % in kg is the rest mass
m_n = 0.26*mass_o;                  % in kg is effective mass of electrons
tempK = 300;                         % in Kelvin (temperature at room)


%program that will model the random motion of the electrons. 
particle_count = 30;      %number of particles
particle_toPlot = 10;      %number of particles to plot in simulation
numItr =100;
Xlength = 2e-7;            %horizontal - based on figure2 in assignment1 given     
Yheight = 1e-7;            %vertical   - based on figure2 in assignment1 given  


%The thermal velocity vth - Assume T = 300 K. 
V_thermal = sqrt(2*k_B*tempK/m_n);       % in m/s 

%spacial step be smaller than 1/100 of the region which is 2e-7, 
% 1/100 *2e-7 = 2e-9 / Vth
dt_step_size = Yheight/(V_thermal*100);    %smallest step allowed
time_step = 500;                       %Simulate for nominally 1000 timesteps.                        
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

%ADD particle positional variable for MFP is given the value of orignal
%position variable
xcord2=xcord;
ycord2=ycord;

%ADD - scatter
dt_scat = dt_step_size/V_thermal;           %new variable for dt - has smae step size as dt_step_size
Pscat_exp = 1 - exp(-dt_scat /t_mn);        %exponial probability scatter equation


%ADD Add in the inner rectangle "bottle neck" boundaries 
%initial ones = use rectangle function to manke BC - based on figure in
%assingment 1 on the dimension of the bottle neck
for i = 1:particle_count
    
    r1=[1 6];
    r2=[1 10];
    % top box defined x-dimensions from 0.8 to 1.2 and y-dimensions from 0.6 to 1
    if xcord(i) > 0.6e-7 && xcord(i) < 1.2e-7 && ycord(i) > 0.6e-7
        
        xcord(i) = randi(r1)/10e-7; %divide th by the wide of the plot region
        ycord(i) = randi(r2)/10e-7;
    end

    % bottom box defined x-dimensions from 0.8 to 1.2 and y-dimensions from 0 to 0.4
    if xcord(i) > 0.6e-7 && xcord(i) < 1.2e-7 && ycord(i) < 0.4e-7
        
        xcord(i) = randi(r1)/10e-7;      %creates an n -by- n codistributed matrix of uniformly distributed random integers
        ycord(i) = randi(r2)/10e-7;
    end
    
    
end
        
       


%before generating values in the loop set all values  to zero
tempMaterial = 0;                     
Vel_bounds_Avg = 0;                    
total_temp = 0;
tempMaterial_prev = 0;
num_Electron_path=0;
distance_between=0;
distance_total=0;
count=0;




%loop through to get calculed value output up to a max threshold
while count < max_time 
    
    x_cordprev = xcord;
    y_cordprev = ycord;
    
    xcord(1:particle_count) = xcord(1:particle_count) + (dt_step_size .* x_vel);  
    ycord(1:particle_count) = ycord(1:particle_count) + (dt_step_size .* y_vel);
    
    
    
    %loop through to check x and y boundries in the region
    for i=1:particle_count
      
        %x boundary
       if xcord(i) > Xlength
           x_cordprev(i) = x_cordprev(i)-Xlength;
           xcord(i) = xcord(i)-Xlength;       
       end
       if xcord(i) < 0
           x_cordprev(i) = x_cordprev(i) + Xlength;
           xcord(i) =xcord(i) +Xlength ;      
       end
       
       %y boundary
       if ycord(i) >Yheight
           y_vel(i) = -y_vel(i);
       end
       if ycord(i)< 0
           y_vel(i) = -y_vel(i);
       end
        
       %ADD in the bottlenect to the trajectory electron traces
       % top box defined x-dimensions from 0.8 to 1.2 and y-dimensions from 0.6 to 1
       %electron collision at the sides of box (top)
       if xcord(i) > 0.8e-7 && xcord(i) < 1.2e-7 && ycord(i) > 0.6e-7
           if ycord(i) > 0.6e-7
               x_vel(i) = -x_vel(i);
           end
           if ycord(i) == 0.6e-7 
               y_vel(i) = -y_vel(i);
           end
       end

     % bottom box defined x-dimensions from 0.8 to 1.2 and y-dimensions from 0 to 0.4
     %electron collision at the sides of box (bottom)
       if xcord(i) > 0.8e-7 && xcord(i) < 1.2e-7 && ycord(i) < 0.4e-7
           if ycord(i) < 0.4e-7 
               x_vel(i) = -x_vel(i);
           end
           if ycord(i) == 0.4e-7 
               y_vel(i) = -y_vel(i);
           end
       end       
       
       
%        %ADD - scatter
%        sdev=sqrt(k_B*tempK/m_n); 
%        dt_scat = dt_step_size/V_thermal;       %new variable for dt - has smae step size as dt_step_size
%        Pscat_exp = 1 - exp(-dt_scat /t_mn);    %exponial probability scatter equation
%     
%        
%        start_scatter = Pscat_exp > randn(particle_count,1); %check condition to start scatter
%        xcord = xcord+x_vel*dt_step_size;
%        ycord = ycord+y_vel*dt_step_size;
%        x_vel(start_scatter)=sdev*randn(sum(start_scatter),1);
%        y_vel(start_scatter)=sdev*randn(sum(start_scatter),1);
       
%        
%         %ADD the MFP   
%         %add the sum of the scatter to find the value of the number of path 
%         num_Electron_path = num_Electron_path + sum(start_scatter); 
%         distance_between = sqrt((xcord2(start_scatter)-xcord(start_scatter)).^2 + ((ycord2(start_scatter)-ycord(start_scatter)).^2));
%         %add the  sum of the distance between particles to the total distance 
%         distance_total = distance_total + sum(distance_between);
%         xcord2(start_scatter) = xcord(start_scatter);      
        
    end   
    
        %calculate the MFP values
        mean_fp_frm_plot = sum(distance_between());  %sum of the distance values in the matrix
        mean_fp_frm_plot_avg = distance_total /num_Electron_path;
        t_mn_frm_plot = mean_fp_frm_plot_avg/V_thermal;
        
        
        %Store valvue of the previos and current positions- for plotting
%       x = [x_cordprev(i) xcord(i)];   
%       y = [y_cordprev(i) ycord(i)]; 

        %do the temperature calc  an check bounds
       Vel_bounds_Avg = mean(sqrt(x_vel.^2 + y_vel.^2));       
       total_temp = (m_n*Vel_bounds_Avg.^2)/(2*k_B);
     
       
       %ADD plot the rectangle bottle nect onto the previos plor
        figure(1) 
        hold on
        subplot(2,1,1)
        axis([0 Xlength 0 Yheight]);
        title(['Q3 Electron enhancement  (' ,num2str(particle_toPlot), ' electrons)']);
        xlabel('xVelocity (nm)');
        ylabel('yVeloncity (nm)');
        rectangle('Position', [0.8e-7 0 0.4e-7 0.4e-7]);
        rectangle('Position', [0.8e-7 0.6e-7 0.4e-7 0.4e-7]);
        hold on;
   
         
        for z=1:particle_toPlot
           %display in subplot 1 the moving electrons             
           colors = hsv(particle_toPlot);          %set to have randoms colors of the stream
           plot([x_cordprev(z) xcord(z)],[y_cordprev(z) ycord(z)], '-', 'markers', 1, 'color', colors(z,:), 'MarkerFaceColor', colors(z,:));
           
        end
       
       %ADD
       
        
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
     
        
        
       pause(0.001); 
       count=count+dt_step_size;        %increment the count one step 
       

end   

%ADD electron density plot - using hist3 function to plot
figure (2);
bin1=30;
bin2=20;
e_v=[xcord',ycord'];
e_density= hist3(e_v, [bin1 bin2]);
title('Q3 - Electron Density Plot ')

% %ADD PLOT the distribution of speeds among the particles at a temperature
% figure (2);
% nbins=25;
% Vel_bounds_Avg= sqrt(x_vel.^2 + y_vel.^2);
% 
% %the distribution of speeds among the particles at a temperature
% histogram(Vel_bounds_Avg,nbins);
% title(' Q2 - Maxwell Boltzman Dist')
% xlabel('Velocity of ELectrons (all directions)(m/s)')
% ylabel('Amount of Electrons ')  