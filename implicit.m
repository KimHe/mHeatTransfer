
YH = 59;% n grids and has n interior points per dimension
XW = 59;
x = linspace(0,0.2,XW+1); %Material Property
y = linspace(0,0.1,YH+1); %Material Property
dx = x(2)-x(1);
h = dx^2/4;
dy =  y(2)-y(1);

p = 1600;% kg/m3,
k = 3;% W/m-K,
c = 0.8;%  J/kg-K.
a = k/(p*c);

Iint = 2:XW+1;
Jint = 2:YH+1;
B = zeros(XW,YH);

Told = zeros(XW+1,YH+1);

Told(1,1:YH) =0;   %BOTTOM
Told(XW,1:YH+1) = -200*x*x'+400*x+300; %TOP;
Told(1:XW+1,1) = 300;  %LEFT
Told(1:XW+1,YH) = 300;  %RIGHT


Tsoln = Told;
B(:,1) = B(:,1)  + Told(Iint,1)*a/h^2; % boundary valuues are added from the left hand side matrix including k and h
B(:,YH) = B(:,YH) + Told(Iint,YH)*a/h^2;
B(1,:) = B(1,:) + Told(1,Jint)*a/h^2;
B(XW,:) = B(XW,:) + Told(XW,Jint)*a/h^2;

% Matrix formations
F = reshape(B,XW*XW,1);
I = speye(XW);
e = ones(XW,1);
T = spdiags([e -4*e e],[-1 0 1],XW,XW);
S = spdiags([e e],[-1 1],XW,XW);
A = (kron(I,T) + kron(S,I))* a / h^2; % A is a tridiagonal spurse matix multiplied and dvided by h^2

Tvec = abs(mldivide(A,F));
Tsoln(Iint,Jint) = reshape(Tvec,XW,XW);
sol = max(max(Tvec)) ;

% Ploting
figure,pcolor(x,y,Tsoln), shading interp

title('Temperature (Implict)'),xlabel('x'),ylabel('y'),colorbar

