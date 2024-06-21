clear; clc;

% Degond's equation: 2D case
% 2 walls boundary condition

% grid
dx = 0.2;
dy = dx;
L = 5;      % [-L,L]^2
nX = floor(2 * L/dx);
nY = floor(2 * L/dy);
x = dx/2-L:dx:L-dx/2;
y = dy/2-L:dy:L-dy/2;

% coefficients
CFL = 0.5;

% Poisson-EQMOM N's
N = 4;

% Degond
tau = 1;
nu = 1;
d = 0.2;
v0 = 1;


RM = zeros(N+1,nY,nX); % 所有结果数据, m0,m1,m2

% Initialize
RM(1,:,:) = 1;
for i = 1:nX
    for j = 1:nY
        theta = (rand() - 0.5) * 2 * pi;    % Random
        eqN = 50;
        ths = (0:1:eqN-1) * 2 * pi / eqN;
        [rmeq_init] = DVM_eq(cos(theta),sin(theta),ths,d,nu);
        rmeq_init = rmeq_init / sum(rmeq_init) * 1;
        for n = 1:N    
            RM(n+1,j,i) = sum( rmeq_init(:) .* exp(1i * n * ths(:)) );
        end
    end
end

disp(sum(RM(1,:,:),'all'));


rho = zeros(N,nY,nX);
phi = zeros(N,nY,nX);
r = zeros(nY,nX);
shift = zeros(nY,nX);
for i = 1:nX
    for j = 1:nY
        M = RM(1:N+1,j,i);
        [rho(:,j,i), phi(:,j,i), r(j,i), shift(j,i) ] = Invn_protect(M);
    end
end


k = 0; tmoment = 0;

F_int_r = zeros(N+3,nY,nX); % m_{-1} to m_{N+1}
F_int_l = zeros(N+3,nY,nX);
F_int_u = zeros(N+3,nY,nX);
F_int_d = zeros(N+3,nY,nX);

flux_X = zeros(N+1,nY,nX+1);% 0 to N
flux_Y = zeros(N+1,nY+1,nX);

RMtp = zeros(N+1,nY,nX);
tsp = 500;

save_k = 0;
result_rho = zeros(nY,nX,1);
result_U = zeros(nY,nX,1);

% iteration starts
while( tmoment<tsp - 1e-6 )
    
    maxRM2 = max(abs(RM(2,:,:)),[],'all');
    km = min(maxRM2 * nu / 2/ d, N);
    dt1 = tau / (maxRM2 * nu * km - d * km^2);
    dt = CFL * min([dx/v0, dy/v0, dt1]);
    
    tmoment = tmoment + dt;
    disp([tmoment,dt,maxRM2]);
    
    
    % flux
    for i = 1:nX
        for j = 1:nY
            for m = -1:N+1
                
                if m == 0
                    F_int_r(m+2,j,i) = - shift(j,i) / 2;
                    F_int_l(m+2,j,i) = - shift(j,i) / 2;
                    F_int_u(m+2,j,i) = - shift(j,i) / 2;
                    F_int_d(m+2,j,i) = - shift(j,i) / 2;
                else
                    F_int_r(m+2,j,i) = - shift(j,i) / 2 / pi * (exp(m * 1i * pi / 2) - exp(-m * 1i * pi / 2))/ 1i / m;
                    F_int_l(m+2,j,i) = - shift(j,i) / 2 / pi * (exp(-m * 1i * pi / 2) - exp(m * 1i * pi / 2))/ 1i / m;
                    F_int_u(m+2,j,i) = - shift(j,i) / 2 / pi * (exp(m * 1i * pi) - 1)/ 1i / m;
                    F_int_d(m+2,j,i) = - shift(j,i) / 2 / pi * (1 - exp(m * 1i * pi))/ 1i / m;
                end

                for n = 1:N
                    [Fnp,Fnn] = Flux(m,phi(n,j,i),r(j,i));
                    F_int_r(m+2,j,i) = F_int_r(m+2,j,i) + rho(n,j,i) * Fnp;
                    F_int_l(m+2,j,i) = F_int_l(m+2,j,i) + rho(n,j,i) * Fnn;
                    [Fnp,Fnn] = Flux(m,phi(n,j,i)-pi/2,r(j,i));
                    F_int_u(m+2,j,i) = F_int_u(m+2,j,i) + rho(n,j,i) * Fnp  * 1i^m;
                    F_int_d(m+2,j,i) = F_int_d(m+2,j,i) + rho(n,j,i) * Fnn  * 1i^m;
                end
                
            end
        end
    end
    
    for i = 2:nX 
        for j = 1:nY
            flux_X(:,j,i) = F_int_l(3:N+3,j,i) + F_int_l(1:N+1,j,i) + F_int_r(3:N+3,j,i-1) + F_int_r(1:N+1,j,i-1);
        end
    end
    for j = 1:nY
        flux_X(:,j,1) = F_int_l(3:N+3,j,1) + F_int_l(1:N+1,j,1) ...
            + ((-1).^(1:1:N+1))' .* conj( F_int_l(3:N+3,j,1) + F_int_l(1:N+1,j,1) );
        flux_X(:,j,nX+1) = F_int_r(3:N+3,j,nX) + F_int_r(1:N+1,j,nX) ...
            + ((-1).^(1:1:N+1))' .* conj( F_int_r(3:N+3,j,nX) + F_int_r(1:N+1,j,nX) );
    end
    
    
    for i = 1:nX
        for j = 2:nY
            flux_Y(:,j,i) = F_int_d(3:N+3,j,i) - F_int_d(1:N+1,j,i) + F_int_u(3:N+3,j-1,i) - F_int_u(1:N+1,j-1,i);
        end
        flux_Y(:,1,i) = F_int_d(3:N+3,1,i) - F_int_d(1:N+1,1,i) + F_int_u(3:N+3,nY,i) - F_int_u(1:N+1,nY,i);
        flux_Y(:,nY+1,i) = F_int_d(3:N+3,1,i) - F_int_d(1:N+1,1,i) + F_int_u(3:N+3,nY,i) - F_int_u(1:N+1,nY,i);
    end
    
    RMtp = RM(1:N+1,:,:) - v0/2*dt/dx*(flux_X(:,:,2:nX+1)-flux_X(:,:,1:nX)) ...
        - v0/2/1i*dt/dy*(flux_Y(:,2:nY+1,:)-flux_Y(:,1:nY,:));
    
    for i = 1:nX
        for j = 1:nY
            % collision
            Omega = RMtp(2,j,i) / abs(RMtp(2,j,i));
            M = RMtp(1:N+1,j,i);
            [rho(:,j,i), phi(:,j,i), r(j,i)] = Invn_protect(M);
            RMtpNp2 = r(j,i)^(N+1) * ( dot(rho(:,j,i), exp(1i * (N+1) * phi(:,j,i))) );
            A = eye(N+1);
            for m = 1:N
                A(m+1,m+1) = 1 + dt/tau*d*m^2;
                A(m+1,m) = -dt/tau*nu*m/2 * Omega;
            end
            for m = 2:N
                A(m,m+1) = dt/tau*nu*(m-1)/2 * conj(Omega);
            end
            
            b = RMtp(:,j,i);
            b(N+1) = b(N+1) - dt/tau*nu*N/2*conj(Omega) * RMtpNp2;
            RM(1:N+1,j,i) = A \ b;                        
            
           
            % moment inversion
            M = RM(1:N+1,j,i);
            [rho(:,j,i), phi(:,j,i), r(j,i), shift(j,i)] = Invn_protect(M);
        end
    end
    
    
    
    figure(1);
    surf(y,x,squeeze(RM(1,:,:)));

    colorbar;
    axis equal;
    
    ylim([-L,L]);
    xlim([-L,L]);
    disp(sum(RM(1,:,:),'all'));
    
    old_save_k = save_k;
    save_k = ceil(tmoment);
    result_rho(:,:,save_k) = RM(1,:,:);
    result_U(:,:,save_k) = RM(2,:,:);
    if old_save_k < save_k
        name = strcat(datestr(datetime,'yyyymmdd'), "_Vicsek2wall_N", num2str(N), "_X", num2str(nX), ".mat");
        save(name);
    end

end

name = strcat(datestr(datetime,'yyyymmdd'), "_Vicsek2wall_N", num2str(N), "_X", num2str(nX), ".mat");
save(name);
