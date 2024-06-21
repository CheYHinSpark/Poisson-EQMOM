clc,clear;

% density variation on a point


% add results of particle models
particle_data = "particle_data_2wall.mat";
load(particle_data);

figure(1);
subplot(4,1,1);
plot(partical_times, avg_rho_wall_l,LineWidth=1);
legend("Particle", Location="northwest");

namelist = ["20240526_Vicsek2wall_N4_X50";
            "20240517_Vicsek2wall_N8_X50";
            "20240522_Vicsek2wall_N16_X50"];
% namelist = ["20240526_Vortex_N4_X50";
%             "20240502_Vortex_N8_X50";
%             "20240507_Vortex_N16_X50"];

% colors = ["-b";
%           "--r";
%           "-.g"];

Pos = [1,1]; Posname = "$\mathbf{x}=(-5,-5)$";
% Pos = [1,25]; Posname = "$\mathbf{x}=(0,-5)$";
for nl = 1:length(namelist)
    name = namelist(nl);
    load(name);

    Max = length(result_rho(1,1,:));
    
    figure(1);
    subplot(4,1,1 + nl);
    plot(1:1:Max,squeeze(result_rho(Pos(1),Pos(2),:)),"LineWidth",1);
    legend(strcat("$N=",num2str(N),"$"),Interpreter="latex",Location="northwest");
    % hold on;
end

for i = 1:4
subplot(4,1,i);
ylim([0 4]);
ylabel("$\rho(\mathbf{x},t)$","Interpreter","latex","FontSize",10);
set(gca,"LineWidth",1);
set(gca,'FontSize',10,'TickLabelInterpreter','latex');
end

subplot(4,1,1);
title(Posname,Interpreter="latex",FontSize=10);
subplot(4,1,2);
xlim([-12,188]);
subplot(4,1,3);
xlim([-2,198]);
subplot(4,1,4);
xlim([-11,189]);
xlabel("$t$","Interpreter","latex","FontSize",10);


set(gcf,'unit','centimeters','position',[0,0,20,16]);


