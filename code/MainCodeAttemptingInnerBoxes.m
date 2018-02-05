% % Assignment1: 100963181 Sonya Stuhec-Leonard

% Effective mass of electron Melectron = 0.26*m0, m0=rest mass
% nominal size of region is 200nm X 100nm
% Electron modeling


%constants
m0 = 9.109e-31; %in kg from source: https://en.wikipedia.org/wiki/Electron
Melectron = 0.26*m0;
k = physconst('Boltzmann'); %Constants in matlab source: https://www.mathworks.com/help/phased/ref/physconst.html?s_tid=gn_loc_drop
T = 300; % temperature in Kalvin

%define thermal velocity source: https://en.wikipedia.org/wiki/Thermal_velocity
v_th = sqrt(k*T/Melectron);

a1 = 3; %acceleration
F1= m0*a1;
numP = 10; %number of particles
iterations = 50; %number of iterations

%box definitions
xmax = 200e-9;
xmin = 0;
ymax = 100e-9;
ymin = 0;

%bottleneak defineitions
TopboxYmax = ymax;
TopboxYmin = 75e-9;

BoxXmax = 125e-9;
BoxXmin = 75e-9;
BboxYmax = 50e-9;
BboxYmin = ymin;

%initalize randome positions for particles
%randmone number between 0 and 100nm or 200nm to be within box
xlocations = rand(numP, 1).*xmax;
ylocations = rand(numP, 1).*ymax;
positions = [xlocations, ylocations];

%generate a randome inital angle for each particle in rad
angle = rand(numP, 1).*2*pi;

% randome velocity angle with magnitude v_th
velocityX = v_th.* cos(angle);
velocityY = v_th.* sin(angle);
velocity = [velocityX, velocityY];

% Randome velocity and gle and magnitude based on Maxell boltzman
% disibution https://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/Kinetics/Rate_Laws/Gas_Phase_Kinetics/Maxwell-Boltzmann_Distributions

% MBfunc = @(v) (Melectron/(2*pi*k*T))^(1/2)*exp((-Melectron*v^2)/(2*k*T));% @(c) (4*pi.*c^2)*(Melectron/(2*pi*k*T))^(3/2)*exp(-Melectron.*c^2/(2*k*T))
% %vels = vecotor of numP ranome velocities.
%
% for index = 1:legnth(vels)
%     weight(index) =  MBfunc(index)
% end
%
% RandVelX = randsample(vels,numP,true,weight)
% RandVelY = randsample(vels,numP,true,weight)

%use 100 steps to get across the region 200nm long
t = (200e-9/v_th)/100;

%inner boxes inital conditions test
% for initx = 1:numP
%     if (positions(initx, 1)<BoxXmax && positions(initx, 1)>BoxXmin) && (positions(initx, 2)> TopboxYmin || positions(initx, 2)< BboxYmax) %for top inner box
%         initx=1;
%         %choose new x particale inital  location
%         positions(initx, 1) = rand(1).*xmax;
%         positions(initx, 2) = rand(1).*ymax;
%         
%     end
% end

for inity = 1:numP
    while (1)
        if positions(inity, 1)<=BoxXmax && positions(inity, 1)>=BoxXmin && (positions(inity, 2)<= TopboxYmin || positions(inity, 2)>= BboxYmax) %for inner box
            %inity=inity-1;
            
            %choose new y particale inital  location
            positions(inity, 1) = rand(1).*xmax;
            %choose new y particale inital  location
            positions(inity, 2) = rand(1).*ymax;
        else
            break
        end
    end
end

for iter =1:1:iterations
    
    %Keep position and velocity form previouse iteration
    oldP = positions;
    oldV = velocity;
    
    %Boundary conditions
    for j=1:2
        for n=1:length(xlocations)
            positions(n, j) = positions(n, j) + velocity(n, j)*t;
            %restrications of x-cordinate of each particle
            if j == 1 %(for all x cordinates)
                if positions(n, 1) <= xmin
                    
                    positions(n, 1) = xmax + velocity(n, 1)*t;
                    
                elseif positions(n, 1)>= xmax
                    
                    positions(n, 1)= xmin + velocity(n, 1)*t;
                    
%                 elseif positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin %within x region of bottle neck
%                     if positions(n, 2) > TopboxYmin %side of top portion in bottle
%                         velocity(n, 1) = -1*velocity(n, 1);
%                     end
%                 elseif positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
%                     if positions(n, 2) == TopboxYmin %top of bottle
%                         velocity(n, 2) = -1*velocity(n, 2);
%                     end
%                 elseif positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
%                     if positions(n, 2)<BboxYmax %side of bottom of bottle
%                         velocity(n, 1) = -1*velocity(n, 1);
%                     end
%                 elseif positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
%                     if positions(n, 2)==BboxYmax
%                         velocity(n, 2) = -1*velocity(n, 2);
%                     end
                    
                end
            end
            
            %y parmaters of region 100X200nm
            if j == 2
                if positions(n, j) <= ymin || positions(n, j) >= ymax
                    velocity(n, 2) = -1*velocity(n, 2);%just negate y component
                    
                elseif positions(n, j) > TopboxYmin %side of top portion in bottle
                    if positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin %within x region of bottle neck
                        velocity(n, 1) = -1*velocity(n, 1);
                    end
                elseif positions(n, j) == TopboxYmin %top of bottle
                    if positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
                        velocity(n, 2) = -1*velocity(n, 2);
                    end
                elseif positions(n, j)<BboxYmax %side of bottom of bottle
                    if positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
                        velocity(n, 1) = -1*velocity(n, 1);
                    end
                elseif positions(n, j)==BboxYmax
                    if positions(n, 1)<=BoxXmax && positions(n, 1)<=BoxXmin
                        velocity(n, 2) = -1*velocity(n, 2);
                    end
                end

            end
           
        end
    end
    %Plot x and y cordinated of position (.b indicates curent posioion
    %only) to plot lines save previouse point and plot the two points as
    %a line. Hold on will keep each segment plotted before
    
    %to add colout put the plot i a loop, make an array of coulour, adn
    %for each plot of particle a select colout n form teh colour
    %array---
    
    %     dim = [.2 .5 .3 .3];
    % Temperature formula from: https://en.wikipedia.org/wiki/Thermal_velocity
    temp = (mean(velocity(:, 1))^2 + mean(velocity(:, 2))^2)*Melectron/k;
    
    Temperature = 'Temperature:';
    
    string = strcat(Temperature, ' ' , num2str(temp));
    %characters = char(Temperature)
    
    figure (1)
    plot(positions(:, 1), positions(:, 2), '.b')
    hold on
    %         h = text(string)
    %         delete(h)
    axis([xmin, xmax, ymin, ymax])
    %     dim1 = [BoxXmin/xmax (BboxYmin+11e-9)/ymax (BoxXmax-BoxXmin)/xmax (BboxYmax-BboxYmin)/ymax];
    %     annotation('rectangle',dim1,'Color','k')
    %     dim2 = [BoxXmin/xmax (TopboxYmin-7.5e-9)/ymax (BoxXmax-BoxXmin)/xmax (TopboxYmax-TopboxYmin)/ymax];
    %     annotation('rectangle',dim2,'Color','k')
    
    pause(0.2)
    title ('Electron Simulation of Trajectories')
    
    
end

%electron denisty map, loop over x, loop over y (starting dimentsions for the density.
% and then loop over particles to see where they are in the region

