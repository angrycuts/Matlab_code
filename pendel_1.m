%the pendulum
ang=0:0.01:2*pi; 
time =0:0.1:pi; 
r=0.04;

%theta = time; 
L = 1;
x0 =0; y0=1; 
w = 0.01; 


figure;  
for i = 1:32
    theta = time(i); %with euler, theta = fun(time(i))
    x = L*cos(-theta) + x0; 
    y = L*sin(-theta) + y0;

    xdata=[x-w, x+w, w, -w];
    ydata=[y, y, L, L]; 
    patch(xdata, ydata, [0,0,1] );
    
    xp=r*cos(ang);
    yp=r*sin(ang);
    patch(x+xp,y+yp, [1,0,0]);
    
    axis([-L L -L L]);
    axis equal;
    hold off; grid on; 
     
    pause(1/1000);
end


%plot(L*cos(ang),L*sin(ang)+L, 'k');  
