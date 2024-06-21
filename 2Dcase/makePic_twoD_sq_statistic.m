clear,clc;

% 6.3 periodic case, error statistic


namelist = ["20240511_Vicsek_N4_X50";
            "20240510_Vicsek_N8_X50";   
            "20240511_Vicsek_N16_X50"];
colors = ["-b";
          "--r";
          "-.g"];

for nl = 1:length(namelist)
    name = namelist(nl);
    load(name);

    % result_U is \rho U
    
    Tn = length(result_U(1,1,:));
    Tn = 501; % 
    Us = zeros(1, Tn - 1);
    Errors_rho = zeros(1, Tn - 1);
    Errors_U = zeros(1, Tn - 1);
    for t = 1:Tn - 1
        Errors_rho(t) = norm(result_rho(:,:,t) - 1) / sqrt(nX * nY);
        Us(t) = sum(result_U(:,:,t),"all") / nX / nY;
        Errors_U(t) = norm(result_U(:,:,t) - Us(t)) / sqrt(nX * nY);
    end
        
    figure(1);
    plot(1:Tn-1, log10(Errors_rho(1:Tn-1)),colors(nl),"LineWidth",2);
    hold on;

    figure(2);
    plot(1:Tn-1, log10(Errors_U(1:Tn-1)),colors(nl),"LineWidth",2);
    hold on;
end

%% fig1, rho
figure(1);
xlabel("$t$","Interpreter","latex","FontSize",20);
ylabel("$\log_{10}\|\rho - \overline{\rho}\|_{L^2_x([-5,5]^2)}$","Interpreter","latex","FontSize",20);
legend(["$N=4$";"$N=8$";"$N=16$"],"Interpreter","latex","Location","southwest");

title("Errors of $\rho$", "Interpreter", "latex", FontSize=20);

set(gca,"LineWidth",2);
set(gca,'FontSize',18,'TickLabelInterpreter','latex');

set(gcf,'unit','centimeters','position',[1,2,17,14]);

%% fig2, U
% rho U
figure(2);
xlabel("$t$","Interpreter","latex","FontSize",20);
ylabel("$\log_{10}\|\rho \mathbf{U} - \overline{\rho \mathbf{U}}\|_{L^2_x([-5,5]^2)}$","Interpreter","latex","FontSize",20);
legend(["$N=4$";"$N=8$";"$N=16$"],"Interpreter","latex","Location","southwest");

title("Errors of $\rho \mathbf{U}$", "Interpreter", "latex", FontSize=20);

set(gca,"LineWidth",2);
set(gca,'FontSize',18,'TickLabelInterpreter','latex');

set(gcf,'unit','centimeters','position',[1,2,17,14]);