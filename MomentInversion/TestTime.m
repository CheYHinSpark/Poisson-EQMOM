clear,clc;

% 202406
% Use this script to test the runtime of different methods. 


Nt = 512;
ths = (0:1:Nt-1) * 2 * pi / Nt;

load("random_data.mat");

N = 32;
t_total = 0;

for i = 1:10000
    R = r(i);
    P = squeeze(p(:,i));
    f = F1(ths, R, P);
    
    % get moments
    Mf = zeros(N+1,1);
    for n = 0:N
        Mf(n + 1) = dot(f, exp(n * 1i * ths)) / Nt * 2 * pi;
    end
    
    % moment inverse
    tic;
    % [rho_i, phi_i, r_i, shift] = Invn(Mf);
    % [rho_i, phi_i, r_i, shift] = Invn_protect(Mf);
    [rho_i, phi_i, r_i, shift] = Invn_noshift(Mf);
    t = toc;
    t_total = t_total + t;
    
end

disp(N);
disp(t_total);



function [f] = F1(theta, R, P)
    f = 0.25 / pi * (1 + ...
        P(1) * sin(theta) + P(2) * sin(2*theta) + P(3) * sin(3 * theta) + ...
        P(4) * sin(5 * theta) + P(5) * sin(7 * theta) + P(6) * sin(11*theta)) + ...
        0.25 / pi * (1 - R^2) ./ (1 - 2 * R * cos(theta) + R^2);
end
