clc,clear;

% 202406
% Compare reuslts with different N's

casename = "case3";

Ns = [8, 12, 16];
taus = [1,0.1,0.01];
colors = ["b--","r-.","k:"];

rhos = zeros(1000,3);
thetas = zeros(1000,3);

runtimes_matrix = zeros(3);

for ns = 1:length(Ns)
    NNe = strcat("N",num2str(Ns(ns)),"e");
for i = 1:length(taus)
    tau = taus(i);
    
    name = strcat(casename,NNe,num2str(tau),"x1000");
    load(strcat(name,".mat"));    
    
    Nx = length(result_rho);
    disp(Nx);
    dx = 10/Nx;
    x = -5+dx/2:dx:5-dx/2;



    figure(1);
    subplot(1,length(taus),i);
    plot(x,result_rho,colors(ns));hold on;


    figure(2);
    subplot(1,length(taus),i);
    plot(x,result_theta,colors(ns));hold on;
    

    
    runtimes_matrix(ns,i) = runtime;
end

    rhos(:,ns) = result_rho(:);
    thetas(:,ns) = result_theta(:);
end

disp(runtimes_matrix);


figure(1);
subplot(1,3,1);
ylabel("$\rho$",Interpreter="latex");
ylim([0.4 1.6]);
xlabel("$x$",Interpreter="latex");
title("$\tau=1$",Interpreter="latex");
subplot(1,3,2);
ylim([0.4 1.6]);
xlabel("$x$",Interpreter="latex");
title("$\tau=0.1$",Interpreter="latex");
subplot(1,3,3);
ylim([0.4 1.6]);
xlabel("$x$",Interpreter="latex");
title("$\tau=0.01$",Interpreter="latex");
legend(["$N=8$","$N=12$","$N=16$"],Interpreter="latex",Location="northwest");

figure(2);
subplot(1,3,1);
ylabel("$\theta$",Interpreter="latex");
ylim([-1.5 1.5]);
xlabel("$x$",Interpreter="latex");
title("$\tau=1$",Interpreter="latex");
subplot(1,3,2);
ylim([-1.5 1.5]);
xlabel("$x$",Interpreter="latex");
title("$\tau=0.1$",Interpreter="latex");
subplot(1,3,3);
ylim([-1.5 1.5]);
xlabel("$x$",Interpreter="latex");
title("$\tau=0.01$",Interpreter="latex");
legend(["$N=8$","$N=12$","$N=16$"],Interpreter="latex",Location="northwest");


disp(norm(rhos(:,1)-rhos(:,2)) * sqrt(dx));
disp(norm(rhos(:,1)-rhos(:,3)) * sqrt(dx));
disp(norm(rhos(:,2)-rhos(:,3)) * sqrt(dx));
