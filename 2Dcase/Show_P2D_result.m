clc,clear;


% 20230703_P2DDFL_N4_X100 用的 invn_protect 0 100

name = "result_P2D_DFL\20230807_P2DDFL_N16_X200.mat";
name = "result_P2D_DFL\20230806_P2DDFL_N8_X200.mat";
% name = "result_P2D_DFL\20230801_P2DDFL_N4_X1000.mat";
% name = "result_P2D_DFL\20230727_P2DDFL_N4_X400.mat";
name = "result_P2D_DFL\20230727_P2DDFL_N4_X200.mat";


% PN 结果
% name = "20240228_P2DDFL_N8_X400_PN.mat";
% name = "20240226_P2DDFL_N8_X200_PN.mat";
% name = "20240226_P2DDFL_N4_X200_PN.mat";
name = "20240225_P2DDFL_N16_X200_PN.mat";
name = "20240312_P2DDFL_N8_X5000_PN.mat";
% name = "20240319_P2DDFL_N8_X2500_PN.mat";
% name = "20240308_P2DDFL_N8_X2000_PN.mat";



% name = "20230808_P2DVicsek_N8_X200.mat";
load(name);                             


Nx = 10 / dx;
x = -5 + dx/2: dx: 5 - dx/2;
y = zeros(1,Nx);

Max = length(result_rho(1,:));


for i = 500:1000
    disp(i);
    figure(1);
    plot(x,result_rho(:,i));hold on; 
    U = real(result_U(:,i)) ./ result_rho(:,i);
    V = imag(result_U(:,i)) ./ result_rho(:,i);
    quiver(x',y',U,V);hold off;
    axis([-5,5,-inf,inf]);
    title(strcat("\it t = ", num2str(tsp / Max * i)));
    pause(0.05);
    xlabel("$x$", "Interpreter","latex");
    ylabel("$\rho$", "Interpreter","latex");
end
% legend(["$\Delta x=0.05$","$\Delta x=0.025$","$\Delta x=0.01$"],"Interpreter","latex");
% legend(["$N=16$","$N=8$","$N=4$"],"Interpreter","latex");
% 
px = 1;
r_s = r(px);
rho_s = rho(:,px);
phi_s = phi(:,px);
shift_s = shift(px);
vMnn = 1000;
rm_N_poisson = zeros(vMnn,1);
ths = (0:1:vMnn-1) * 2 * pi / vMnn;
for i=1:vMnn
    for k = 1:N
        rm_N_poisson(i) = rm_N_poisson(i) + (1-r_s^2)/(2*pi) * rho_s(k) / (1 + r_s^2 - 2*r_s * cos(phi_s(k) - ths(i)));
    end
    rm_N_poisson(i) = rm_N_poisson(i) - shift_s / 2 / pi;
end
figure(2);
plot(ths,rm_N_poisson);

