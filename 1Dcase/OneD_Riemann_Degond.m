clear; clc;

% 202406
% Degond's equation: 1-D Riemann problem

% grid
dx = 0.01;
L = 5;      % [-L, L]
nnod = floor(2 * L/dx);
x = dx/2-L:dx:L-dx/2;

% coefficients
CFL = 0.5;

% Degond's equation
tau = 1;
nu = 1;
d = 0.2;
v0 = 1;

% Poisson-EQMOM N's
N = 12;

% data
tsp = 4;
RM = zeros(N+2,nnod);   % m0,m1,...,m_N,m_N+1

% initianlize
ndisc0 = floor(nnod/2);

% 3 settings in REF: Gamba Motsch 2014
% case 1
% casename = "case1";
% r1 = 2;
% th1 = 1.7;
% r2 = 0.218;
% th2 = 0.5;

% case 2
% casename = "case2";
% r1 = 1;
% th1 = 1.5;
% r2 = 2;
% th2 = 1.83;

% case 3
casename = "case3";
r1 = 1;
th1 = 1;
r2 = 1;
th2 = -1;

eqN=200;
ths = (0:1:eqN-1) * 2 * pi / eqN;

[rmeq_init] = DVM_eq(cos(th1),sin(th1),ths,d,nu);
rmeq_init = rmeq_init / sum(rmeq_init) * r1;
RM(1,1:ndisc0) = r1;
for i = 1:N
    RM(i+1,1:ndisc0) = sum( rmeq_init(:) .* exp(1i * i * ths(:)) );
end

[rmeq_init] = DVM_eq(cos(th2),sin(th2),ths,d,nu);
rmeq_init = rmeq_init / sum(rmeq_init) * r2;
RM(1,ndisc0+1:nnod) = r2;
for i = 1:N
    RM(i+1,ndisc0+1:nnod) = sum( rmeq_init(:) .* exp(1i * i * ths(:)) );
end

rho = zeros(N,nnod);
phi = zeros(N,nnod);
r = zeros(nnod,1);
shift = zeros(nnod,1);
for i = 1:nnod
    M = RM(1:N+1,i);
    [rho(:,i), phi(:,i), r(i), shift(i) ] = Invn(M);
    RM(N+2,i) = r(i)^(N+1) * ( dot(rho(:,i), exp(1i * (N+1) * phi(:,i))) );
end

% start iteration
k = 0; tmoment = 0;
Fp = zeros(N+3,nnod); Fn = Fp; % m_{-1} to m_N+1
flux = zeros(N+1,nnod+1);
RMtp = zeros(N+1,nnod);

tic;
while( tmoment < tsp - 1e-6 )
    k = k+1;
    
    km = min(nu / 2/ d, N);
    dt1 = tau / (nu * km - d * km^2);
    dt = CFL * min([dx/v0, dt1]);
    
    tmoment = tmoment + dt;
    disp([tmoment,dt]);

    for i = 1:nnod
       
       for m = -1:N+1
           
           if m == 0
               Fp(m+2,i) = - shift(i) / 2;
               Fn(m+2,i) = - shift(i) / 2;
           else
               Fp(m+2,i) = - shift(i) / 2 / pi * (exp(m * 1i * pi / 2) - exp(-m * 1i * pi / 2))/ 1i / m;
               Fn(m+2,i) = - shift(i) / 2 / pi * (exp(-m * 1i * pi / 2) - exp(m * 1i * pi / 2))/ 1i / m;
           end
           
           for n = 1:N
               [Fnp,Fnn] = Flux(m,phi(n,i),r(i));
               Fp(m+2,i) = Fp(m+2,i) + rho(n,i) * Fnp;
               Fn(m+2,i) = Fn(m+2,i) + rho(n,i) * Fnn;
           end           

       end
    end
    
    for i = 2:nnod 
       flux(:,i) = Fp(3:N+3,i-1) + Fp(1:N+1,i-1) + Fn(3:N+3,i) + Fn(1:N+1,i);
    end
    flux(:,1) = Fp(3:N+3,1) + Fp(1:N+1,1) + Fn(3:N+3,1) + Fn(1:N+1,1);
    flux(:,nnod+1) = Fp(3:N+3,nnod) + Fp(1:N+1,nnod) + Fn(3:N+3,nnod) + Fn(1:N+1,nnod);
    
    RMtp = RM(1:N+1,:)-v0/2*dt/dx*(flux(:,2:nnod+1)-flux(:,1:nnod));
    

    for i = 1:nnod
        % collision term
        Omega = RMtp(2,i) / abs(RMtp(2,i));
        M = RMtp(1:N+1,i);
        [rho(:,i), phi(:,i), r(i)] = Invn(M);
        RMtpNp2 = r(i)^(N+1) * ( dot(rho(:,i), exp(1i * (N+1) * phi(:,i))) );
        A = eye(N+1);
        for m = 1:N
            A(m+1,m+1) = 1 + dt/tau*d*m^2;
            A(m+1,m) = -dt/tau*nu*m/2 * Omega;
        end
        for m = 2:N
            A(m,m+1) = dt/tau*nu*(m-1)/2 * conj(Omega);
        end
        b = RMtp(:,i);
        b(N+1) = b(N+1) - dt/tau*nu*N/2*conj(Omega) * RMtpNp2;
        RM(1:N+1,i) = A \ b;
    end
    
    
    % moment inversion
    for i = 1:nnod
        M = RM(1:N+1,i);
        [rho(:,i), phi(:,i), r(i), shift(i)] = Invn(M);
        RM(N+2,i) = r(i)^(N+1) * ( dot(rho(:,i), exp(1i * (N+1) * phi(:,i))) );
    end
end
runtime = toc;

result_rho = RM(1,:);
plot(x,RM(1,:));

result_theta = angle(RM(2,:));

save(strcat(casename,"N",num2str(N),"e",num2str(tau),"x",num2str(nnod),".mat"),...
    "N","tau","dx","result_rho","result_theta","d","nu","runtime","k");







