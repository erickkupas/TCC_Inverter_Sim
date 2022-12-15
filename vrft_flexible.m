%Flexible VRFT Function
%Chrystian Lenon Remes - 2016, Aug, 25th
%Based on Ricardo Scheid Filho

function [rho,eta,varargout] = vrft_flexible(uk1,uk2,yk1,yk2,Td0,Cbar,rho0,enpd2,itermax,tol)

%uk is a vector, containg experiment input data
%yk is a vector, containg experiment output data
%Td is the Desired Transfer Function, defined by the user (must be proper)
%nz if the number of zeros to be find for T(z,n)
%C_bar is the controller structure, defined by the user
%p_o (rho_o) is the initial parameter rho to be defined by the user. Must result in a stable loop
%Ta is the Sampling Period
%iter_max is the number of iterations to be executed in order to find rho and eta
%Cclass: informs if the desired controller is at the class

Td = tf(Td0);
Ta = Td.Ts;
numTd = Td.num{1};
denTd = Td.den{1};  %Td Denominator
nz = length(roots(numTd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Verification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(uk1) ~= length(yk1)) %Verify if uk and yk have the same lengths
    error('Input u(k) and Output y(k) Vectors must have the same length')
end

if (length(uk2) ~= length(yk2)) %Verify if uk and yk have the same lengths
    error('Input u(k) and Output y(k) Vectors must have the same length')
end

if (length(uk1) ~= length(uk2)) %Verify if uk and yk have the same lengths
    error('Inputs u1(k) and u2(k) Vectors must have the same length')
end

if (isproper(Td0) == 0) %Verify if Td is proper
    error('Td(z) must be a proper Transfer Function (deg(num) < deg(den))')
end

if (nz < 0) %Verify if the number of zeros is positive
    error('The number of zeros to be found must be a non-negative number')
end

if (isproper(Cbar) == 0) %Verify if C_bar is proper
    error('C_bar(z) must be a proper Transfer Function (deg(num) < deg(den))')
end

if (length(Cbar) ~= length(rho0)) %Verify if C_bar and p_o have the same length
    error('C_bar(z) and p_o must have the same length')
end

if (itermax < 1) %Verify if iter_max >= 1
    error('The number of iterations must be greater or equal to 1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%verify if the number of zeros is OK: T(z,n) is proper.
if ((length(denTd)-1) < nz)
    error('The chosen number of zeros "nz" of T(z,n) gives a non-proper TF. Choose nz <= deg(denTd)')
end

%Verify if identification of pd2 of Td is enabled
if enpd2 == 1
    lambda = roots(numTd);
    pd1 = max(abs(roots(denTd)));
    pd2 = (lambda*(1-pd1))/(lambda-pd1);
    if (pd2 > pd1^4)
        pd2 = pd1^4;
    end
    Td = tf(zpk(lambda,[pd1 pd2],(1-pd1)*(1-pd2)/(1-lambda),Ta));
    denTd = Td.den{1};
end

%Create F(z) = [F1(z);F2(z);...Fn(z)] vector
F = tf(1,denTd,Ta);     %F(1) = F1
if (nz+1)>= 2           %if the number of zeros is >= 1,
    for ii = 2:(nz+1)   %create a vector of transfer functions F(z)
        F = [tf([1 zeros(1,ii-1)],denTd,Ta) ; F]; %F(2..n) = [F2 ... Fn]
    end
end

%Compute the initial controller from p_o and defined controller structure C_bar
Cz_vrft = minreal(tf(rho0'*Cbar),1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VRFT Flexible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th_ant = [rho0; numTd(2:end)'];
dth = 1;
iter = 1

while (iter <= itermax) && (dth > tol)
    
    Lz = minreal(Td*(1-Td));
    phi_wn = lsim(F*Lz,uk1) + lsim(F*Lz*Cz_vrft,yk1);
    phi_wn_iv = lsim(F*Lz,uk2) + lsim(F*Lz*Cz_vrft,yk2);
    phi_y2n = lsim(Lz*Cz_vrft,yk1);
    
    eta = ((phi_wn_iv'*phi_wn)^-1)*phi_wn_iv'*phi_y2n;    %Find ni (eta) that minimizes the cost function (Least Squares)
    lambda = -eta(2)/eta(1);
    
    if (enpd2 == 1) && (lambda > 1)
        pd2 = (lambda*(1-pd1))/(lambda-pd1);
        if (abs(pd2) > pd1^4)
            pd2 = pd1^4;
        end
        Td = tf(zpk(lambda,[pd1 pd2],(1-pd1)*(1-pd2)/(1-lambda),Ta));
        denTd = Td.den{1};
        for ii=1:length(F)
            F(ii).den{1} = denTd;
        end
    else
        Td = sum(denTd)/(1-lambda)*tf([1 -lambda],denTd,Ta);
    end
    
    Lz = minreal((eye(1) - Td));    %Find the new L(z) filter from T(z,ni)
    phi_p = lsim(Cbar*Lz*(eye(1)-Td),yk1);    %Apply C_bar*L(z)*(1-T(z,ni)) into yk
    zeta_p_iv = lsim(Cbar*Lz*(eye(1)-Td),yk2);
    phi_u = lsim(Td*Lz,uk1);       %Apply T(z,ni)*L(z) into yk
    
    rho = ((zeta_p_iv'*phi_p)^-1)*zeta_p_iv'*phi_u; %Find pi+1 (rho) that minimizes the cost function (Least Squares)
    
    th = [rho; eta];
    dth = (th-th_ant)'*(th-th_ant)/(length(th));
    th_ant = th;
    
    Cz_vrft = tf(minreal(rho'*Cbar));   %Find the Controller Cz_vrft for the next iteration
    iter = iter+1
end

varargout{1} = dth;
if enpd2 == 1
    varargout{2} = pd2;
end

end