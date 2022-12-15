%Standard VRFT Function
%Chrystian Lenon Remes - 2016, Aug, 25th

function rho = vrft_standard(uk,yk1,yk2,Td,C_bar)

%uk is a vector, containg experiment input data
%yk is a vector, containg experiment output data
%Td is the Desired Transfer Function, defined by the user (must be proper)
%C_bar is the controller structure, defined by the user
%Ta is the Sampling Period

if (length(uk) ~= length(yk1)) %Verify if uk and yk have the same lengths
    error('Input u(k) and Output y(k) Vectors must have the same length')
end

if (length(uk) ~= length(yk2)) %Verify if uk and yk have the same lengths
    error('Input u(k) and Output y(k) Vectors must have the same length')
end

if (isproper(Td) == 0) %Verify if Td is proper
    error('Td(z) must be a proper Transfer Function (deg(num) < deg(den))')
end

if (isproper(C_bar) == 0) %Verify if C_bar is proper
    error('C_bar(z) must be a proper Transfer Function (deg(num) < deg(den))')
end

%Compute the filter L(z)
Lz = (eye(1)-Td);     %From VRFT Background

uLk = lsim(Td*Lz,uk);      %Apply L(z) filter into uk signal
phiL = lsim(Lz*C_bar*(eye(1)-Td),yk1);  %apply L(z) filter into obtained from yk
zetaL = lsim(Lz*C_bar*(eye(1)-Td),yk2);  %apply L(z) filter into obtained from yk
    
rho = ((zetaL'*phiL)^-1)*zetaL'*uLk; %Find the parameters p (rho) of the controller through Least Squares