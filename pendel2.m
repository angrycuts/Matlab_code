%% Clear everything
clear all;
close all;
clc

%% Declaring variables

%WindowSize
wSize=[-2 2 -0.6 2];
wWidth = 0.1;


%Materials
radius = [0.1 0.051];
volume = [4*pi*(radius(1)^3)/3 4*pi*(radius(2)^3)/3];
area = [4*pi*radius(1)^2 4*pi*radius(2)^2];
density = [11340 11340];
mass = [density(1)*volume(1) density(2)*volume(2)];

%Pendulum variables
rope_length = [0.5 0.5];
rope_width = 0.01; 
circle_vec = 0:0.01:2*pi;

%Constants
air_constant = [0.47 0.47];
g = 9.82;
tStep = 0.01;

% Initial angle = pi/4 & initial velocity & acceleration= 0
theta_zero = [-(pi/2) (pi/2)];
w_zero = [30 -33];
a_zero = [0 0];

%% Pendulum
theta = [theta_zero(1) theta_zero(2)];
velocity = [w_zero(1) w_zero(2)];
acceleration = [a_zero a_zero];

figure;
%Calculating and drawing the pendulum
for i=1:100
    clf;
    
    %Calc air-resistance
    airres = [(0.5*(velocity(1)^2)*area(1)*air_constant(1))/mass(1) (0.5*(velocity(2)^2)*area(2)*air_constant(2))/mass(2)];
    
    %Direction of air-resistance depends on velvector
    if(velocity(1) < 0)
        airres(1) = airres(1)*-1;
    end
    if(velocity(2) < 0)
        airres(2) = airres(2)*-1;
    end
    
    %Calc Acc/Vel/Pos
    acceleration = [-(g/rope_length(1))*sin(theta(1)), -(g/rope_length(2))*sin(theta(2))];
    velocity = [euler(tStep, velocity(1), acceleration(1))-airres(1), euler(tStep, velocity(2), acceleration(2))- airres(2)];
    theta = [euler(tStep, theta(1), velocity(1)) euler(tStep, theta(2), velocity(2))];
    
    %Position of the spehere
    x = [rope_length(1)*sin(theta(1))-1 1+rope_length(2)*sin(theta(2))];  
    y = [rope_length(1)*(1-cos(theta(1))) rope_length(2)*(1-cos(theta(2)))];
    
    %Drawing the rope
    xdata=[x(1)-rope_width, x(1)+rope_width, rope_width-1, -rope_width-1; x(2)-rope_width, x(2)+rope_width, 1+rope_width, 1-rope_width];
    ydata=[y(1), y(1), rope_length(1), rope_length(1); y(2), y(2), rope_length(2), rope_length(2)]; 
    patch(xdata(1,:), ydata(1,:), [0,0,1]);
    patch(xdata(2,:), ydata(2,:), [0,0,1]);
    
    %Drawing the spehere
    xp = radius(1)*cos(circle_vec);
    xp2 = radius(2)*cos(circle_vec);
    yp = radius(1)*sin(circle_vec);
    yp2 = radius(2)*sin(circle_vec);
    patch(x(1)+xp,y(1)+yp, [1,0,0]);
    patch(x(2)+xp2,y(2)+yp2, [0,1,0]);
    
    %Drawing the Wall
    wallXData=[wSize(1), wSize(1)+wWidth, wSize(1)+wWidth, wSize(2)-wWidth, wSize(2)-wWidth, wSize(2), wSize(2), wSize(1)];
    wallYData=[wSize(4), wSize(4), wSize(3)+wWidth,wSize(3)+wWidth, wSize(4),wSize(4), wSize(3), wSize(3)]; 
    patch(wallXData, wallYData, [1,0,0] );
    
    %WindowSize
    xlim([wSize(1) wSize(2)]);
    ylim([wSize(3) wSize(4)]);

    grid on; 
    
    pause(1/100);
end

% Projectile

% Initial values
xAcc = [0 0];
yAcc = [acceleration(1)*sin(theta(1))-g, acceleration(2)*sin(theta(2)) - g];
xVel = [velocity(1)*cos(theta(1)) velocity(2)*cos(theta(2))];
yVel = [velocity(1)*sin(theta(1)) velocity(2)*sin(theta(2))];
xPos = [x(1) x(2)];
yPos = [y(1) y(2)];

%calculating and drawing projectile
for i = 1:700
    clf;
    
    %Air resistance;
    xAirres = [(0.5*(xVel(1)^2)*area(1)*air_constant(1))/mass(1),  (0.5*(xVel(2)^2)*area(2)*air_constant(2))/mass(2)];
    yAirres = [(0.5*(yVel(1)^2)*area(1)*air_constant(1))/mass(1),  (0.5*(yVel(2)^2)*area(2)*air_constant(2))/mass(2)];      
    %Calculate Velocity
    xVel = [euler(tStep, xVel(1), xAcc(1))- xAirres(1), euler(tStep, xVel(2), xAcc(2)) - xAirres(2)]; 
    yVel = [euler(tStep, yVel(1), yAcc(1))- yAirres(1), euler(tStep, yVel(2), yAcc(2)) - yAirres(2)];
    %Calculate Posistion
    xPos = [euler(tStep, xPos(1), xVel(1)), euler(tStep, xPos(2), xVel(2))]; 
    yPos = [euler(tStep, yPos(1), yVel(1)), euler(tStep, yPos(2), yVel(2))];
    
    %Collide    
    if(sqrt((xPos(1) - xPos(2))^2 + (yPos(1) - yPos(2))^2) < (radius(1)+radius(2))) 
        r2 = radius(1)+radius(2);
        
        % Set VelVec to the old ones.
        yVel = tempyVel;
        xVel = tempxVel;
        
        % Calc ang
        angleOne = atan2(yVel(1),xVel(1));
        angleTwo = atan2(yVel(2),xVel(2));
        
        %Move the objects from eachother
        move = 0.001;
        while(sqrt((xPos(1) - xPos(2))^2 + (yPos(1) - yPos(2))^2) < r2)
        xPos(1) = xPos(1) - (move * cos(angleOne));
        yPos(1) = yPos(1) - (move * sin(angleOne));
        
        xPos(2) = xPos(2) - (move * cos(angleTwo));
        yPos(2) = yPos(2) - (move * sin(angleTwo));
        end
        
        %Find nor.plane
        x = [xPos(2) yPos(2)] - [xPos(1) yPos(1)];
        x = x/norm(x);
        
        v1 = [xVel(1) yVel(1)];
        x1 = dot(x, v1); 
        v1x = x*x1;
        v1y = v1 - v1x;
        
        x = -x;
        v2 = [xVel(2) yVel(2)];
        x2 = dot(x, v2);
        v2x = x*x2;
        v2y = v2-v2x;
        
        totMass = mass(1) + mass(2);
        
        newV1 = (v1x * ((mass(1)-mass(2)) /  totMass)) + (v2x*  ((mass(2)) /totMass))     + v1y;
        newV2 = (v1x * ((mass(1))       /  totMass)) + (v2x*  ((mass(2)-mass(1))/totMass))+ v2y;
        
        xVel = [newV1(1) newV2(1)];
        yVel = [newV1(2) newV2(2)];
    end
    
    %Hit Ground
    if(yPos(1) < -0.5+radius(1))
        yVel(1) = -yVel(1)*0.8;
        yPos(1) = -0.5+radius(1);
        xVel(1) = xVel(1)*0.8;
    end
    
    if(yPos(2) < -0.5+radius(2))
        yVel(2) = -yVel(2)*0.8;
        yPos(2) = -0.5+radius(2);
        xVel(2) = xVel(2)*0.8;
    end
    
    %Hit Wall
    if(xPos(1)-radius(1) < wSize(1)+wWidth || xPos(1)+radius(1) > wSize(2)-wWidth)
        xVel(1) = -xVel(1);
    end
    if(xPos(2)-radius(2) < wSize(1)+wWidth || xPos(2)+radius(2) > wSize(2)-wWidth)
        xVel(2) = -xVel(2);
    end
     
    %Drawing the spehere
    xp=radius(1)*cos(circle_vec);
    xp2 = radius(2)*cos(circle_vec);
    yp=radius(1)*sin(circle_vec);
    yp2=radius(2)*sin(circle_vec);
    
    patch(xPos(1)+xp,yPos(1)+yp, [1,0,0]);
    patch(xPos(2)+xp2,yPos(2)+yp2, [0,1,0]);

    %Drawing the wall
    wallXData=[wSize(1), wSize(1)+wWidth, wSize(1)+wWidth, wSize(2)-wWidth, wSize(2)-wWidth, wSize(2), wSize(2), wSize(1)];
    wallYData=[wSize(4), wSize(4), wSize(3)+wWidth,wSize(3)+wWidth, wSize(4),wSize(4), wSize(3), wSize(3)]; 
    patch(wallXData, wallYData, [0,0,0] );
    
    %WindowSize
    xlim([wSize(1) wSize(2)]);
    ylim([wSize(3) wSize(4)]);
    grid on; 
    
    tempxPos = xPos;
    tempyPos = yPos;
    tempyVel = yVel;
    tempxVel = xVel;
    
    
    pause(1/100);
end