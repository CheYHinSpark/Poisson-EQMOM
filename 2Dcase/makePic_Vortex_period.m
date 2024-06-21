clc,clear;

% 该脚本是为了统计4面墙算例中某一个点上的密度变化情况
 
% namelist = ["20240526_Vicsek2wall_N4_X50";
%             "20240517_Vicsek2wall_N8_X50";
%             "20240522_Vicsek2wall_N16_X50"];
namelist = ["20240526_Vortex_N4_X50";
            "20240502_Vortex_N8_X50";
            "20240507_Vortex_N16_X50"];

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
    subplot(3,1,nl);
    plot(1:1:Max,squeeze(result_rho(Pos(1),Pos(2),:)),"LineWidth",1);
    legend(strcat("$N=",num2str(N),"$"),Interpreter="latex",Location="northwest");
    % hold on;
end

for i = 1:3
subplot(3,1,i);
ylim([0 5]);
ylabel("$\rho(\mathbf{x},t)$","Interpreter","latex","FontSize",10);
set(gca,"LineWidth",1);
set(gca,'FontSize',10,'TickLabelInterpreter','latex');
end

subplot(3,1,1);
title(Posname,Interpreter="latex",FontSize=10);
xlim([0,200]);
yticks([0,1,2,3,4,5]);
subplot(3,1,2);
xlim([0,200]);
yticks([0,1,2,3,4,5]);
subplot(3,1,3);
xlim([0,200]);
yticks([0,1,2,3,4,5]);
xlabel("$t$","Interpreter","latex","FontSize",10);


set(gcf,'unit','centimeters','position',[0,0,20,13]);


