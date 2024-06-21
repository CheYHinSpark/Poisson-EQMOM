function [rho_eq] = DVM_eq(xO,yO,ths,d,nu)

    % get the von Mises distribution. 
    % (x0,y0):  the velocity
    % ths:      the theta's discrete vector
    % d, nu:    coefficients

    N = length(ths);
    d_th = 2 * pi / N;
    xs = cos(ths);
    ys = sin(ths);
    
    mO = sqrt(xO^2 + yO^2);
    xO = xO / mO;
    yO = yO / mO;

    QU = zeros(N);

    for i = 1:N
        QU(i,i) = -2 * d / d_th^2;
    end

    for i = 1:N-1
        QU(i,i+1) = d / d_th^2 + nu * (xO * ys(i+1) - yO * xs(i+1)) / 2 / d_th;
    end
    QU(N,1) = d / d_th^2 + nu * (xO * ys(1) - yO * xs(1)) / 2 / d_th;

    for i = 1:N-1
        QU(i+1,i) = d / d_th^2 - nu * (xO * ys(i) - yO * xs(i)) / 2 / d_th;
    end
    QU(1,N) = d / d_th^2 - nu * (xO * ys(N) - yO * xs(N)) / 2 / d_th;

    rho_eq = null(QU);

    sum_rho = sum(rho_eq);

    rho_eq = rho_eq / sum_rho;

end
