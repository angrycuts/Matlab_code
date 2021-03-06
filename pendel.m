%% Clear everything
clear all;
close all;
clc

%% Declaring variables

%WindowSize
wSize=[-2 2 -0.6 2];
wWidth = 0.1;

%Pendulum variables
radius = 0.04;
volume = 4*pi*(radius^3)/3;
area = 4*pi*radius^2;
density = 11340;
mass = density*volume;
rope_length = 0.5;
rope_width = 0.01; 
circle_vec = 0:0.01:2*pi;

%Constants
air_constant = 0.47;
g = 9.82;
tStep = 0.01;

% Initial angle = pi/4 & initial velocity & acceleration= 0
theta_zero = -(pi/2);
w_zero = 70;
a_zero = 0;

%Triangle
sSTri = zeros(2,51);
sSTri(1,:) = [1:0.01:1.5];
sSTri(2,:) = 0.5;

%% Pendulum
theta = theta_zero;
velocity = w_zero;
acceleration = a_zero;

figure;
%Calculating and drawing the pendulum
for i=1:100
    clf;
    airres = (0.5*(velocity^2)*area*air_constant)/mass;
    if(velocity < 0)
        airres = airres*-1
    end
    
    acceleration = -(g/rope_length)*sin(theta)
    velocity = euler(tStep, velocity, acceleration)- airres;
    theta = euler(tStep, theta, velocity);
    
    %Position of the spehere
    x = rope_length*sin(theta);  
    y = rope_length*(1-cos(theta));
    
    %Drawing the rope
    xdata=[x-rope_width, x+rope_width, rope_width, -rope_width];
    ydata=[y, y, rope_length, rope_length]; 
    patch(xdata, ydata, [0,0,1] );
    
    %Drawing the spehere
    xp=radius*cos(circle_vec);
    yp=radius*sin(circle_vec);
    patch(x+xp,y+yp, [1,0,0]);
    
    %Drawing the triangle
    triXData = [1.25, 1.5, 1];
    triYData = [1, 0.5, 0.5];
    patch(triXData, triYData, [0,1,0]);
    
    %Drawing the Wall
    wallXData=[wSize(1), wSize(1)+wWidth, wSize(1)+wWidth, wSize(2)-wWidth, wSize(2)-wWidth, wSize(2), wSize(2), wSize(1)];
    wallYData=[wSize(4), wSize(4), wSize(3)+wWidth,wSize(3)+wWidth, wSize(4),wSize(4), wSize(3), wSize(3)]; 
    patch(wallXData, wallYData, [1,0,0] );
    
    %WindowSize
    xlim([wSize(1) wSize(2)]);
    ylim([wSize(3) wSize(4)]);
    grid on; 
    
    pause(1/1000);
end

% Projectile

xAcc = 0;
yAcc = acceleration*sin(theta) - g;
xVel = velocity*cos(theta);
yVel = velocity*sin(theta);
xPos = x;
yPos = y;

%calculating and drawing projectile
for i = 1:400
    clf;
    
    xAirres = (0.5*(xVel^2)*area*air_constant)/mass;
    yAirres = (0.5*(yVel^2)*area*air_constant)/mass;
    
    xVel = euler(tStep, xVel, xAcc) - xAirres; 
    yVel = euler(tStep, yVel, yAcc) - yAirres;
    
    xPos = euler(tStep, xPos, xVel);
    yPos = euler(tStep, yPos, yVel);
    
    %Calculating the bouncing angle
    angle = atan(yVel/xVel);
    
    %Hit the triangle
    if((xPos < 1.5)  && (xPos > 1.0) && (yPos > 0.5-radius) && yPos < 0.5)
        
    end

    %Hit Ground
    if(yPos < -0.5+radius)
        yVel = -yVel*0.8;
        yPos = -0.5+radius;
        xVel = xVel*0.8;
    end
    
    %Hit Wall
    if(xPos < wSize(1)+wWidth || xPos > wSize(2)-wWidth)
        xVel = -xVel;
    end
    
    %Drawing the spehere
    xp=radius*cos(circle_vec);
    yp=radius*sin(circle_vec);
    patch(xPos+xp,yPos+yp, [1,0,0]);
    
    %Drawing the triangle
    triXData = [1.25, 1.5, 1];
    triYData = [1, 0.5, 0.5];
    patch(triXData, triYData, [0,1,1]);

    %Drawing the wall
    wallXData=[wSize(1), wSize(1)+wWidth, wSize(1)+wWidth, wSize(2)-wWidth, wSize(2)-wWidth, wSize(2), wSize(2), wSize(1)];
    wallYData=[wSize(4), wSize(4), wSize(3)+wWidth,wSize(3)+wWidth, wSize(4),wSize(4), wSize(3), wSize(3)]; 
    patch(wallXData, wallYData, [0,0,0] );
    
    %WindowSize
    xlim([wSize(1) wSize(2)]);
    ylim([wSize(3) wSize(4)]);
    grid on; 
    
    pause(1/1000);
end