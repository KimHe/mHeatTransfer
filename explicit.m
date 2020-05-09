
clear %Generating Grid
YH = 60;
XW = 60;
x = linspace(0,0.2,XW); %Material Property
dx = x(2)-x(1);
y = linspace(0,0.1,YH); %Material Property
dy = y(2)-y(1); ; %1/YH;



p = 1600;% kg/m3,
k = 3;% W/m-K,
c = 0.8;%  J/kg-K.
a = k/(p*c);


T = zeros(XW);
T(1,1:YH) = 0;   %BOTTOM
T(XW,1:YH) = 200*x*x'+400*x+300; %TOP;
T(1:XW,1) = 300;  %LEFT
T(1:XW,YH) = 300;  %RIGHT
dt = dx^2/4;
TOL = 1e-3;
error = 1; k = 0;
    while error > TOL
        k = k+1;
          Told = T;
          for i = 2:XW-1
              for j = 2:YH-1
              T(i,j) =Told(i,j)+a*dt*((Told(i+1,j)-2*Told(i,j)+Told(i-1,j))/dx^2 ...
                      + (Told(i,j+1)-2*Told(i,j)+Told(i,j-1))/dy^2);
              end
          end
          T(XW,1:YH) = -200*x*x'+400*x+300; %TOP;
          error = max(max(abs(Told-T)));
    end

figure,pcolor(x,y,T),shading interp,

title('Temperature (explict)'),xlabel('x'),ylabel('y'),colorbar

