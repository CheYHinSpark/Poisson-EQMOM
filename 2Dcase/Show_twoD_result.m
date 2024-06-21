clc,clear;

% name = "result_vortex/20230613_vortexN6.mat";
% name = "result_periodic/20230622_DFL_N8_X50.mat";
% name = "result_periodic/20230607_xxxxDFL4_sq2.mat";   % 这个看起来最好
% name = "xxxxDFL4_sq2_1.5pi.mat";
name = "20240517_Vicsek2wall_N8_X50.mat";

load(name);

                                                                                             
Nx = 2 * L / dx;
x = -L + dx/2: dx: L - dx/2;
y = -L + dx/2: dx: L - dx/2;

% surf(y,x,squeeze(RM(1,:,:)));
Max = length(result_rho(1,1,:));

i = 1;
while (i <= Max)
    disp(i);
    figure(1);
    surf(y,x,result_rho(:,:,i));hold on;
    % colormap(hot);
    u=real(result_U(:,:,i))./result_rho(:,:,i);
    v=imag(result_U(:,:,i))./result_rho(:,:,i);
    quiver(x,y,u,v);hold off;
    axis equal;
    axis([-L,L,-L,L,0,5]);
    title(strcat("\it t = ", num2str(tsp / Max * i)));
    pause(0.1);
    i = i + 1;
end

% 
% px = 7;
% py = 3;
% r_s = r(px,py);
% rho_s = rho(:,px,py);
% phi_s = phi(:,px,py);
% shift_s = shift(px,py);
% vMnn = 1000;
% rm_N_poisson = zeros(vMnn,1);
% ths = (0:1:vMnn-1) * 2 * pi / vMnn;
% for i=1:vMnn
%     for k = 1:N
%         rm_N_poisson(i) = rm_N_poisson(i) + (1-r_s^2)/(2*pi) * rho_s(k) / (1 + r_s^2 - 2*r_s * cos(phi_s(k) - ths(i)));
%     end
%     rm_N_poisson(i) = rm_N_poisson(i) - shift_s / 2 / pi;
% end
% 
% plot(ths,rm_N_poisson);


