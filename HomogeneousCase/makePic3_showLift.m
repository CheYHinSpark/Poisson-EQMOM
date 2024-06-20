clc,clear;

% Show the effect of lifting

N = 8;
d = 0.8;
nu = 1;

NN = 1024;

ths = (0:1:NN-1) * 2 * pi / NN;

lifts = [0,0.5];
line_color = ["r--","b-","k:","m-."];
legends = ["";"";"";"\it M"];


theta = pi / 4;
[f_init] = DVM_eq(cos(theta),sin(theta),ths,d,nu);

RM = zeros(N+1,1);
RM(1) = sum(f_init);
for i = 1:N
    RM(i+1) = dot( f_init(:) , exp(1i * i * ths(:)) );
end

poisson = zeros(NN,1);
    
for j = 1:length(lifts)
    lift = lifts(j);
    M = RM;
    M(1) = M(1) + lift;
    [ rho, phi, r, shift ] = Invn_noshift( M );
    
    for i = 1:NN
        poisson(i,j) = - lift / 2 / pi;
        for k = 1:length(rho)
            fpk = rho(k) * 1 / 2 / pi * (1 - r^2) / (1 - 2 * r .* cos(ths(i) - phi(k)) + r^2);
            poisson(i,j) = poisson(i,j) + fpk;
        end
    end
    
    plot(ths, poisson(:,j), line_color(j));
    hold on;
    
    legends(j) = strcat("\it l = ", num2str(lift));
end

% Invn.m
[ rho, phi, r, shift ] = Invn( RM );
for i = 1:NN
    poisson(i,j+1) = - shift / 2 / pi;
    for k = 1:length(rho)
        fpk = rho(k) * 1 / 2 / pi * (1 - r^2) / (1 - 2 * r .* cos(ths(i) - phi(k)) + r^2);
        poisson(i,j+1) = poisson(i,j+1) + fpk;
    end
end
plot(ths, poisson(:,j+1), line_color(j+1));
hold on;
legends(j+1) = strcat("\it l = ", num2str(shift));


vonMises = f_init * NN / 2 / pi;
plot(ths,vonMises,line_color(j+2));hold on;

xlabel("\it \theta");
ylabel("\it f(\theta)");
xlim([0,2*pi]);
legend(legends);
