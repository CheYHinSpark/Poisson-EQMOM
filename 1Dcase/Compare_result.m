clc,clear;

% 202406
% Compare results with different tau's

casename = "case3";

taus = [1,0.1,0.01];
colors = ["b--","r-.","k:"];

rhos = zeros(1000,3);
thetas = zeros(1000,3);
N=12;

for i = 1:length(taus)
    tau = taus(i);

    name = strcat(casename,"N",num2str(N),"e",num2str(tau),"x1000");
    load(strcat(name,".mat"));    

    Nx = length(result_rho);
    disp(Nx);
    dx = 10/Nx;
    x = -5+dx/2:dx:5-dx/2;



    figure(1);
    plot(x,result_rho,colors(i));hold on;


    figure(2);
    plot(x,result_theta,colors(i));hold on;

    rhos(:,i) = result_rho(:);
    thetas(:,i) = result_theta(:);

end

figure(1);
a_rho = readmatrix(strcat(casename,"\rho_Macro.csv"));
ax = a_rho(:,1);
plot(ax, a_rho(:,2));hold on;
xlim([-5,5]);
xlabel("\it x");
ylabel("\it \rho");
legend(["\it \tau = 1","\it \tau = 0.1","\it \tau = 0.01","Macro"]);


figure(2);
a_theta = readmatrix(strcat(casename,"\theta_Macro.csv"));
ax = a_theta(:,1);
plot(ax, a_theta(:,2));hold on;
xlim([-5,5]);
xlabel("\it x");
ylabel("\it \theta");
legend(["\it \tau = 1","\it \tau = 0.1","\it \tau = 0.01","Macro"]);
