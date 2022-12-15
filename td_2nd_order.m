function [Td,kt,z1,w] = td_2nd_order(p1,p2,Ts,f0)

w = 2*pi*f0; 
z = exp(i*w*Ts); 

syms k z1 real 

A = (z-p1)*(z-p2);
B = k*(z-z1); 

% Solving for k and z1

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
kt = double(k); 

z1 = sol.z1;
z1 = double(z1); 

Td = zpk([z1],[p1 p2],kt,Ts); % Reference Model - VRFT