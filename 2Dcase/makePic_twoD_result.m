clc,clear;

% draw heat maps

% colormap
colormap1 = ones(31,3);
colormap1(1:10,3) = 1.0:-0.1:0.1;
colormap1(11:31,3) = 0;
colormap1(11:20,2) = 1.0:-0.1:0.1;
colormap1(21:31,2) = 0;
colormap1(21:30,1) = 1.0:-0.1:0.1;
colormap1(31,1) = 0;

colormap2 = zeros(21,3);
colormap2(1:21,2) = 1.0:-0.05:0;
colormap2(1:21,1) = 1.0;

% data names
name = "20240511_Vicsek_N16_X50.mat";
name = "20240522_Vicsek2wall_N16_X50.mat";
name = "20240526_Vortex_N4_X50.mat";
load(name);

                                                                                             
Nx = 10 / dx;
x = -5 + dx/2: dx: 5 - dx/2;
y = -5 + dx/2: dx: 5 - dx/2;

Max = length(result_rho(1,1,:));


i = 201;
tmoment = tsp * i / Max;
disp(i);
figure(1);
imagesc(y,x,result_rho(:,:,i));

colormap(colormap1);
% caxis([0.99 1.01]);%clim([0.99 1.01]);
% colormap(colormap2);
cb = colorbar;
hold on

% result_U = squeeze(RM(2,:,:));
U = result_U(:,:,i);
% too many arrows
ddd = 2;
Ux = Nx / ddd;
ux = -5 + dx * ddd / 2 : ddd * dx : 5 - dx * ddd / 2;
u=real(U(:,:)) ./ result_rho(:,:,i);
v=imag(U(:,:)) ./ result_rho(:,:,i);
qu = zeros(Ux,Ux);
qv = zeros(Ux,Ux);
for i = 1:Ux
    for j = 1:Ux
        qu(j,i) = sum(u(ddd*j-ddd+1:ddd*j,ddd*i-ddd+1:ddd*i),"all");
        qv(j,i) = sum(v(ddd*j-ddd+1:ddd*j,ddd*i-ddd+1:ddd*i),"all");
    end
end
quiver(ux,ux,qu,qv,0.6,"LineWidth",2);hold off;
axis([-5,5,-5,5,0,1]);
axis equal;

ax = gca;
ax.YDir = 'normal';
set(cb,"LineWidth",2,'FontSize',18,'TickLabelInterpreter','latex');
set(gca,"LineWidth",2);
set(gca,'FontSize',18,'TickLabelInterpreter','latex');
xlabel("$x$",'Interpreter','latex','FontSize',20);
ylabel("$y$",'Interpreter','latex','FontSize',20);
% title(strcat("$t = ", num2str(floor(tmoment)),"$"),'Interpreter','latex','FontSize',20);
xlim([-5 5]);
ylim([-5 5]);

set(gcf,'unit','centimeters','position',[1,2,16,13.5])
