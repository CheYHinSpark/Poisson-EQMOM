clear; clc;

% Homogeneous problem for Degond's equation

% coefficients
CFL = 0.5;
tau = 1;

% Different N numbers
Ns = [4,8,16,32];

lineStyle = ["-.";"--";":";"-"];

L2err_lists = zeros(1,length(Ns));

M_err_lists = zeros(1,length(Ns));

C0err_lists = zeros(1,length(Ns));

time_cost = zeros(1,length(Ns));


eqN = 1024;
f_last = zeros(eqN,length(Ns));


% for itt = 1:100 % repeat for 100 times to statistic runtime
for n = 1:length(Ns)
    N = Ns(n);
    
    % coefficients for Degond's equation
    nu = 1;
    d = 0.2;
    v0 = 1;
    
    
    tsp = 20;   % time
    t_list = zeros(1);
    RM = zeros(N+2, 1); % m0,m1,...,m_N,m_N+1
    
    % Initialize
    ths = (0:1:eqN-1) * 2 * pi / eqN;
    f_init = (1 + cos(4 * ths)) .* exp( -cos(pi * (ths / 2 / pi + (ths / 2 / pi) .^ 4)));
    rho_init = sum(f_init) * 2 * pi / eqN;
    RM(1) = rho_init;
    for i = 1:N
        RM(i+1) = dot( f_init(:) , exp(1i * i * ths(:)) ) * 2 * pi / eqN;
    end
    

    rho = zeros(N,1);
    phi = zeros(N,1);
    r = 0;
    shift = 0;
    M = RM(1:N+1);
    [rho(:), phi(:), r, shift ] = Invn(M);
    RM(N+2) = r^(N+1) * ( dot(rho(:), exp(1i * (N+1) * phi(:))) );
    
    
    % Iteration starts
    k = 0; tmoment = 0;
    RMtp = zeros(N+1,1);
    while( tmoment < tsp - 1e-8)
        tic;
    
        k = k+1;
    
        km = min(nu / 2/ d, N);
        dt = CFL * tau / (nu * km - d * km^2);
        tmoment = tmoment + dt;
        t_list(k) = tmoment;
        disp([tmoment, dt]);
    
        RMtp = RM(1:N+1,1);
    
        % Collision term
        Omega = RMtp(2,1) / abs(RMtp(2,1));
        theta = angle(Omega);
        RMtpNp2 = r^(N+1) * ( dot(rho(:,1), exp(1i * (N+1) * phi(:,1))) );
        A = eye(N+1);
        for m = 1:N
            A(m+1,m+1) = 1 + dt/tau*d*m^2;
            A(m+1,m) = -dt/tau*nu*m/2 * Omega;
        end
        for m = 2:N
            A(m,m+1) = dt/tau*nu*(m-1)/2 * conj(Omega);
        end
        b = RMtp(:,1);
        b(N+1) = b(N+1) - dt/tau*nu*N/2*conj(Omega) * RMtpNp2;
        RM(1:N+1,1) = A \ b;
    
        % moment inversion
        M = RM(1:N+1,1);
        [rho(:,1), phi(:,1), r, shift] = Invn(M);
        RM(N+2,1) = r^(N+1) * ( dot(rho(:,1), exp(1i * (N+1) * phi(:,1))) );
    
    
        tictoc = toc;
        time_cost(n) = time_cost(n) + tictoc;
    
        % error statistic
        f = zeros(eqN,1);
        for i = 1:eqN
            for j = 1:N
                f(i) = f(i) + rho(j) * (1 - r^2) / (1 + r^2 - 2 * r * cos(phi(j) - ths(i)));
            end
            f(i) = (f(i) - shift) / 2 / pi;
        end
        [VMeq] = DVM_eq(cos(theta),sin(theta),ths,d,nu);
        VMeq = VMeq / sum(VMeq) * rho_init * eqN / 2 / pi;
        RMeq = zeros(N + 2,1);
        for ri = 0:N+1
            RMeq(ri + 1) = dot(VMeq, exp(ri * 1i * ths)) / eqN * 2 * pi;
        end

        L2err_lists(k,n) = norm(f-VMeq) * sqrt(2 * pi / eqN);
        C0err_lists(k,n) = max(abs(f-VMeq),[],"all");
        M_err_lists(k,n) = norm(RM - RMeq) / norm(RMeq);
    end


    %% turn off figures when collecting runtime
    figure(n);
    plot(ths(32:32:eqN), f(32:32:eqN), 'o',LineWidth=2); hold on;
    plot(ths, VMeq, '-',LineWidth=2); hold on;
    plot(ths, f_init, ':',LineWidth=2); hold off;
    xlim([0, 2 * pi]);
    xticks([0,0.5*pi,pi,1.5*pi,2*pi]);
    xticklabels({'0','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2\pi$'});
    ax = gca;
    ax.TickLabelInterpreter = "latex";
    ax.LineWidth = 2;
    ax.FontSize = 20;
    xlabel("$\theta$",Interpreter="latex",FontSize=20);
    ylabel("$f(\theta)$",Interpreter="latex",FontSize=20);
    legend(["$f(t=20,\theta)$";"$M_{\bar{\theta}}(t=20,\theta)$";"$f_0(\theta)$"], Interpreter="latex", FontSize=18,Location="northwest");
    title(strcat('$t=', num2str(tmoment), '$, $N=', num2str(N),'$'), Interpreter="latex", FontSize=20);
    set(gcf,'Units','centimeters','Position',[2 2 15 12]);

    figure(length(Ns) + 1);
    semilogy(t_list,L2err_lists(:,n),lineStyle(n),LineWidth=2); hold on;
    % figure(length(Ns) + 2);
    % plot(t_list,log10(C0err_lists(:,n)),lineStyle(n),LineWidth=2); hold on;
    figure(length(Ns) + 3);
    semilogy(t_list,M_err_lists(:,n),lineStyle(n),LineWidth=2); hold on;
    
end
% end % repeat for 100 times to statistic runtime

%% turn off figures when collecting runtime
figure(length(Ns) + 1);
ylabel("$\|f(\theta) - M_{\bar{\theta}}(\theta)\|_{L^2_{\theta}}$",Interpreter="latex",FontSize=20);
ylim([1e-5,1e1]);
yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1]);
xlabel("$t$",Interpreter="latex",FontSize=20);
legend(["$N=4$";"$N=8$";"$N=16$";"$N=32$"], Interpreter="latex", FontSize=18,Location="southwest");
ax = gca;
ax.TickLabelInterpreter = "latex";
ax.LineWidth = 2;
ax.FontSize = 20;

set(gcf,'Units','centimeters','Position',[2 2 15 12]);

figure(length(Ns) + 3);
ylabel("$|\mathbf{m}[f] - \mathbf{m}[M_{\bar{\theta}}]|/|\mathbf{m}[M_{\bar{\theta}}]|$",Interpreter="latex",FontSize=20);
ylim([1e-6,1e0]);
yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]);
xlabel("$t$",Interpreter="latex",FontSize=20);
legend(["$N=4$";"$N=8$";"$N=16$";"$N=32$"], Interpreter="latex", FontSize=18,Location="southwest");
ax = gca;
ax.TickLabelInterpreter = "latex";
ax.LineWidth = 2;
ax.FontSize = 20;
set(gcf,'Units','centimeters','Position',[2 2 15 12]);