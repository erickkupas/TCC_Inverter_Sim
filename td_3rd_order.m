function [Td,kt,z1,w] = td_3rd_order(p1,p2,p3,z2,Ts,f0)

w = 2*pi*f0; % Freq. angular da rede [rad/s]

z = exp(i*w*Ts); % z na freq. wTs

syms k z1 real % Variáveis simbólicas

% Dividindo em duas equações A e B 

A = (z-p1)*(z-p2)*(z-p3);

B = k*(z-z1)*(z-z2); 

% Separando a parte real da imaginária (2 eq. 2 incógnitas)

reA = real(A);
imA = imag(A);

reB = real(B);
imB = imag(B);

reEq = reA == reB;
reEq(k,z1) = reEq;

imEq = imA == imB;

reEq = subs(reEq(k,z1));

sol = vpasolve([imEq reEq], [k z1]);

k = sol.k;

kt = double(k); % Ganho da Td

z1 = sol.z1;

z1 = double(z1); % Zero 1 da Td

Td = zpk([z1],[p1 p2],kt,Ts); % Modelo de referência Td para o VRFT