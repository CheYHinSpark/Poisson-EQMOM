clc,clear;

% Show the errors with different N's and lifts. 

Ns = [4,8,12,16,24,32];
lifts = [0,0.5,1,2];
d = 0.4;
nu = 1;

NN = 2048;

ths = (0:1:NN-1) * 2 * pi / NN;

line_color = ["r--","b-","k:","m-."];
legends = [];


theta = pi / 4;
[f_init] = DVM_eq(cos(theta),sin(theta),ths,d,nu);
vonMises = f_init * NN / 2 / pi;

poisson = zeros(NN,length(Ns));

L2errs = zeros(length(Ns),length(lifts));
C0errs = zeros(length(Ns),length(lifts));
MNp1errs = zeros(length(Ns),length(lifts));
MNp1_relative_errs = zeros(length(Ns),length(lifts));


for m = 1:length(lifts)
    lift = lifts(m);
    for j = 1:length(Ns)

        N = Ns(j);
        RM = zeros(N+1,1);
        RM(1) = sum(f_init) + lift;
        for i = 1:N
            RM(i+1) = dot( f_init(:) , exp(1i * i * ths(:)) );
        end
        RMNp1 = dot( f_init(:) , exp(1i * (N + 1) * ths(:)) );
        
        [ rho, phi, r ] = Invn_noshift( RM );

        for i = 1:NN
            poisson(i,j) = - lift / 2 / pi;
            for k = 1:length(rho)
                fpk = rho(k) * 1 / 2 / pi * (1 - r^2) / (1 - 2 * r .* cos(ths(i) - phi(k)) + r^2);
                poisson(i,j) = poisson(i,j) + fpk;
            end
        end
        RMNp1_poisson = 0;
        for k = 1:length(rho)
            RMNp1_poisson = RMNp1_poisson + rho(k) * r^(N+1) * exp(1i * (N+1) * phi(k));
        end

        L2err = norm(vonMises - poisson(:,j)) * sqrt(2*pi / NN);
        C0err = max(abs(vonMises - poisson(:,j)));
        L2errs(j,m) = L2err;
        C0errs(j,m) = C0err;
        MNp1errs(j,m) = abs(RMNp1 - RMNp1_poisson);
        MNp1_relative_errs(j,m) = abs(RMNp1 - RMNp1_poisson) / abs(RMNp1);
        

        disp(MNp1_relative_errs(j,m));
    end
    figure(1);
    plot(Ns, log10(L2errs(:,m))); hold on;
end


xlabel("N");
ylabel("log_{10} err");

legend(["lift=0","lift=0.5","lift=1","lift=2"],Location="southwest");
