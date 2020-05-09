function main


    fprintf('\n Please choose the discretization scheme \n');
    fprintf(' ... input "explicit" when you need the explicit scheme \n');
    fprintf(' ... input "implicit" when you need the implicit scheme \n');

    scheme = input('\n Enter ... : ', 's');

    obj = jsonDecode(fileread('config.json'));
    alpha = obj.k / (obj.rho * obj.cp);

    if strcmp(scheme, 'explicit')

        epsilon = 1; k = 0;

        T = zeros(obj.xGrid);

        x = linspace(0, obj.xLength, obj.xGrid);
        y = linspace(0, obj.yLength, obj.yGrid);

        dx = x(2) - x(1);
        dy = y(2) - y(1);
        dt = dx^2 / 4;

        T(1:obj.xGrid, 1) = obj.leftBoundTemp;
        T(1:obj.xGrid, obj.yGrid) = obj.rightBoundTemp;
        T(1, 1:obj.yGrid) = obj.bottomBoundTemp;
        T(obj.xGrid, 1:obj.yGrid) = 200 *x*x' + 400*x + 300; %TOP;


        while epsilon > obj.TOL

            k = k + 1;
            Told = T;

            for i = 2:obj.xGrid-1
                for j = 2:obj.yGrid-1
                    T(i,j) = Told(i,j) + alpha * dt ...
                        * ( (Told(i+1, j) - 2 * Told(i, j) + Told(i-1, j)) / dx^2 ...
                        + ( Told(i, j+1) - 2 * Told(i, j) + Told(i, j-1) ) / dy^2 );
                end
            end

            T(obj.xGrid, 1:obj.yGrid) = -200 * x*x' + 400*x + 300; %TOP;

            epsilon = max(max(abs(Told - T)));

            if mod(k, 5000) == 0
                figure(01), pcolor(x,y,T), shading interp,
            end

        end

        figure(01), pcolor(x,y,T), shading interp,

        title('Temperature (explict)'), xlabel('x length'), ylabel('y length'), colorbar

    elseif strcmp(scheme, 'implicit')


        obj.xGrid = obj.xGrid - 1;
        obj.yGrid = obj.yGrid - 1;

        x = linspace(0, obj.xLength, obj.xGrid);
        y = linspace(0, obj.yLength, obj.yGrid);

        dx = x(2) - x(1);
        dy = y(2) - y(1);
        h = dx^2 / 4;

        Iint = 2:obj.xGrid+1;
        Jint = 2:obj.yGrid+1;
        B = zeros(obj.xGrid, obj.yGrid);

        Told = zeros(obj.xGrid+1, obj.yGrid+1);

        Told(1, 1:obj.yGrid) = obj.bottomBoundTemp;
        Told(obj.xGrid, 1:obj.yGrid+1) = -200 * x*x' + 400*x + 300; %TOP;
        Told(1:obj.xGrid+1, 1) = obj.leftBoundTemp;
        Told(1:obj.xGrid+1, obj.yGrid) = obj.rightBoundTemp;


        Tsoln = Told;
        B(:, 1) = B(:, 1) + Told(Iint, 1) * alpha / h^2; % boundary valuues are added from the left hand side matrix including k and h
        B(:, obj.yGrid) = B(:, obj.yGrid) + Told(Iint, obj.yGrid) * alpha / h^2;
        B(1, :) = B(1, :) + Told(1, Jint) * alpha / h^2;
        B(obj.xGrid, :) = B(obj.xGrid, :) + Told(obj.xGrid, Jint) * alpha / h^2;

        % Matrix formations
        F = reshape(B, obj.xGrid*obj.xGrid, 1);
        I = speye(obj.xGrid);
        e = ones(obj.xGrid, 1);
        T = spdiags([e -4*e e], [-1 0 1], obj.xGrid, obj.xGrid);
        S = spdiags([e e], [-1 1], obj.xGrid, obj.xGrid);
        A = (kron(I, T) + kron(S, I)) * alpha / h^2; % A is a tridiagonal spurse matix multiplied and dvided by h^2

        Tvec = abs(mldivide(A, F));
        Tsoln(Iint, Jint) = reshape(Tvec, obj.xGrid, obj.xGrid);
        sol = max(max(Tvec));

        % Ploting
        figure, pcolor(x,y,Tsoln), shading interp

        title('Temperature (Implict)'), xlabel('x length'), ylabel('y length'), colorbar

    else

        fprintf('ERROR: Please check you have input the corrrect option \n');

    end

end
