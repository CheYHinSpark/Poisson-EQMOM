clc,clear;

% show Poisson kernel

NN = 10000;

ths = (0:1:NN) * 2 * pi / NN;
theta = pi / 3;

r = [0.2,0.4,0.6,0.8];
line_color = ["r--","b-","k:","m-."];
legends = ["";"";"";""];

poisson = zeros(NN+1,1);

for j = 1:length(r)
    rj = r(j);
    for i = 1:NN+1
        poisson(i) = 1 / 2 / pi * (1 - rj^2) / (1 - 2 * rj .* cos(ths(i) - theta) + rj^2);
    end
    
    plot(ths, poisson, line_color(j));
    hold on;
    
    legends(j) = strcat("\it r = ", num2str(rj));
end

xlabel("\it \theta");
ylabel("\it P_r(\theta - \pi/3)");
xlim([0,2*pi]);
legend(legends);
