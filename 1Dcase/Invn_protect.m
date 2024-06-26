function [ rho, phi, r, lift ] = Invn_protect( M )

% 202406
% This inversion algorithm has a limitation on lift. 
% Thus it is more stable than "Invn.m". 
% Use this one when "Invn.m" fails. 

%% Step 1: solve r
N = length(M)-1;

rho = zeros(N,1);
phi = zeros(N,1);

Mr = @(r) M./ power(r,0:N)';
lmin = @(r) real( eigs( toeplitz(Mr(r)), 1, 'smallestreal' ) );

rmax = 1;
r = min(abs(M(2))/M(1), rmax);
if (r == rmax)
    disp(abs(M(2))/M(1));
end

a = r; b = 1;
lift =  -lmin(r);
smax = 100; % limitation on the lift
smin = 0;   % limitation on the lift
while (lift > smax || lift < smin)
    if (lift > smax)
        a = r;
    end
    if (lift < smin)
        b = r;
    end
    r = (a + b) / 2;
    lift = -lmin(r);
end
M(1) = M(1) + lift;


%% Step 2: solve nodes (\phi)
P = zeros(N,N+1);   % P(:,i) are orthonomal polynomials from P_0 to P_{N-1}��
U = zeros(N,N);     % The matrix J
Pzk = zeros(N,N+1); % Pzk(i,j) is the inner product of P_{i-1} and z^{j-1}
P(1,1) = 1 / sqrt(M(1));
Pzk(1,:) = P(1,1) * M(:) ./ (r.^(0:1:N))';
q = zeros(1,N+1);
for k = 0:N-2
    q(1) = 0; q(2:N+1) = P(k+1,1:N);
    for m = 0:k
        U(k+1,m+1) = dot(conj(Pzk(m+1,:)),q);
        q = q-U(k+1,m+1)*P(m+1,:);
    end
    
    for n = 0:k+1
        U(k+1,k+2) = U(k+1,k+2) + conj(q(n+1)) * M(k+2-n) * r^(n-k-1) * q(k+2);
    end
    U(k+1,k+2) = sqrt(U(k+1,k+2));
    P(k+2,:) = q/U(k+1,k+2);
    
    for m = k+1:N
        for n = 0:k+1
            Pzk(k+2,m+1) = Pzk(k+2,m+1) + conj(P(k+2,n+1)) * M(m-n+1) * r^(n-m);
        end
    end
end
q(1) = 0; q(2:N+1) = P(N,1:N);
for m = 0:N-1
    U(N,m+1) = dot(conj(Pzk(m+1,:)),q);
end
Z = eig(U);
phi(1:N) = angle(Z);

%% Step 3: solve weights (w)
Y = zeros(1,N);
for k = 1:N
   Y(k) = M(N+1-k)*r^(k-N);
end
V = vander(exp(1i*phi));
rho = abs(Y/V);
end
