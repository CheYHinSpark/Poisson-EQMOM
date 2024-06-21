clc,clear;

% density variation on a point

% namelist = ["20240526_Vicsek2wall_N4_X50";
%             "20240517_Vicsek2wall_N8_X50";
%             "20240522_Vicsek2wall_N16_X50"];
namelist = ["20240526_Vortex_N4_X50";
            "20240502_Vortex_N8_X50";
            "20240507_Vortex_N16_X50"];

colors = ["-b";
          "--r";
          "-.g"];

Pos = [1,1]; Posname = "$\mathbf{x}=(-5,-5)$";
% Pos = [1,25]; Posname = "$\mathbf{x}=(0,-5)$";
for nl = 1:length(namelist)
    name = namelist(nl);
    load(name);

    Max = length(result_rho(1,1,:));
    
    figure(1);
    plot(1:1:Max,squeeze(result_rho(Pos(1),Pos(2),:)),colors(nl),"LineWidth",2);
    hold on;
end


xlim([0 200]);
xlabel("$t$","Interpreter","latex","FontSize",20);
ylabel("$\rho(\mathbf{x},t)$","Interpreter","latex","FontSize",20);
legend(["$N=4$";"$N=8$";"$N=16$"],"Interpreter","latex","Location","northwest");
title(Posname,"Interpreter","latex","FontSize",20);

set(gca,"LineWidth",2);
set(gca,'FontSize',18,'TickLabelInterpreter','latex');

set(gcf,'unit','centimeters','position',[1,2,18,12]);


