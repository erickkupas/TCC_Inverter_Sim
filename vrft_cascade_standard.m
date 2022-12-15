%VRFT Function for Cascade Loops
%Author: Chrystian L. Remes
%-----------------------------------------------------
%Cascade Loop 1 has Ci(z) at the feedback inner branch
%Cascade Loop 1 has Ci(z) at the direct inner branch
%-----------------------------------------------------
%VRFT_cascade_MP(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,looptype,varargin)
%-----------------------------------------------------
%u1 is the input, yi1 is the inner output, ye2 is the outer output
%u2, yi2 and ye2 are the instrumental variables (IV). For cascade loop 2, u2 is
%not required
%-----------------------------------------------------
%Cibar and Cebar are the controller classes for inner and outer loops
%-----------------------------------------------------
%Tde is the reference model for outer loop
%-----------------------------------------------------
%looptype = 1 selects cascade loop 1, looptype = 2 selects cascade loop 2
%-----------------------------------------------------
%Optional Variables
%enSiz = 0 or 1 -> selects estimation of Si(z) for L(z) filter in cascade 1
%u2 is the IV required for cascade 1
%itermax is the number of iterations for cascade 1. Default itermax=5
%tol is the tolerance for the norm of parameter variation in cascade 1. Default tol=1e-6
%Tdi is the inner reference model required for cascade 2
%Examples
%Cascade 1: vrft_cascade_standard(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,looptype,enSiz,u2,itermax,tol)
%Cascade 1: vrft_cascade_standard(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,looptype,Tdi)

%%
function [rho,varargout] = vrft_cascade_standard(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,looptype,varargin)

%% Consistance verification
if (size(yi1) ~= size(yi2))
    error('yi1 and yi2 must be the same length');
end
if (size(ye1) ~= size(ye2))
    error('ye1 and ye2 must be the same length');
end
if (size(u1) ~= size(yi1))
    error('u and yi must be the same length');
end
if (size(u1) ~= size(ye1))
    error('u and ye must be the same length');
end
if (isproper(Tde) ~= 1)
    error('Tde must be proper');
end
if (isproper(Cebar) ~= 1)
    error('Cebar must be proper');
end
if (isproper(Cibar) ~= 1)
    error('Cibar must be proper');
end
if (nargin < 9)
    error('Not enough inputs for any cascade loop');
end

%% Switch between the cascade loops and verify consistence
switch looptype
    case 1
        switch length(varargin)
            case 1
                enSiz = varargin{1};
                if (enSiz == 0)
                    u2 = u1;
                    itermax = 5;
                    tol = 1e-6;
                    [rho,varargout{1}] = VRFT_Cascade_MP1(u1,u2,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,enSiz,itermax,tol);
                else
                    error('Not enough inputs for Cascade 1 with Si(z) Estimation: make enSiz=0 or provide a Instrumental Variable u2 with enSiz=1');
                end
            case 2
                enSiz = varargin{1};
                u2 = varargin{2};
                itermax = 5;
                tol = 1e-6;
                [rho,varargout{1}] = VRFT_Cascade_MP1(u1,u2,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,enSiz,itermax,tol);
            case 4
                enSiz = varargin{1};
                u2 = varargin{2};
                itermax = varargin{3};
                tol = varargin{4};
                [rho,varargout{1}] = VRFT_Cascade_MP1(u1,u2,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,enSiz,itermax,tol);
            otherwise
                error('VRFT with Cascade 1: Incorrect number of parameters');
        end
    case 2
        if nargin == 10
            Tdi = varargin{1};
            rho = VRFT_Cascade_MP2(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tdi,Tde);
            varargout{1} = [];
        else
            error('VRFT with Cascade 2: Incorrect number of parameters');
        end
    otherwise
        error('Select a loop type: "1" for Cascade Loop 1 (Feedback Ci) and "2" for Cascade Loop 2 (Direct Ci)');
end

end

%% VRFT - Cascade 1
function [rho,Ndth] = VRFT_Cascade_MP1(u1,u2,yi1,yi2,ye1,ye2,Cibar,Cebar,Tde,enSiz,itermax,tol)

Ta = Tde.Ts;

if (size(u1) ~= size(u2))
    error('u1 and u2 must be the same length');
end
if (itermax < 0)
    error('maximum number of iterations must be positive');
end
if (tol < 0)
    error('tolerance must be positive');
end

%first estimation of rho with Si(z) = 1
Siz = tf(1,1,Ta);
L = Siz*(1-Tde);
uL = lsim(L*Tde,u1);
phiL1 = [lsim(L*Cebar*(1-Tde),ye1)'; -lsim(L*Cibar*Tde,yi1)']';
phiL2 = [lsim(L*Cebar*(1-Tde),ye2)'; -lsim(L*Cibar*Tde,yi2)']';
rho1 = (phiL2'*phiL1)^-1*(phiL2'*uL);
if (sum(rho1<0)>=1)
    rho1 = lsqlin(phiL1,uL,[],[],[],[],1e-12*ones(1,length(rho1))',[inf*ones(1,length(rho1))]'); %lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
end

ii = 1;
dth = ones(length(rho1),1)';
Ndth = [];

if (enSiz == 1)
    p_e = length(Cebar);
    while (ii <= itermax) && (norm(dth) > tol)
        Ciz1 = rho1(p_e+1:end)'*Cibar;
        ue1 = lsim(Ciz1,yi1)+u1;
        ue2 = lsim(Ciz1,yi2)+u2;
        data = iddata(ue1,u1);
        Si_id = ivx(data,[2 3 0],ue2);
        if (max(abs(roots(Si_id.B)))>1)
            Si_id = arx(data,[2 3 0]);
            if (max(abs(roots(Si_id.B)))>1)
                Si_id.B = 1;
                Si_id.A = 1;
            end
        end
        Siz = sum(Si_id.B(1))/sum(Si_id.A(1))*tf(Si_id.A,Si_id.B,Ta);

        L = Siz*(1-Tde);
        uL = lsim(L*Tde,u1);
        phiL1 = [lsim(L*Cebar*(1-Tde),ye1)'; -lsim(L*Cibar*Tde,yi1)']';
        phiL2 = [lsim(L*Cebar*(1-Tde),ye2)'; -lsim(L*Cibar*Tde,yi2)']';
        rho1 = (phiL2'*phiL1)^-1*(phiL2'*uL);
        if (sum(rho1<0)>=1)
            rho1 = lsqlin(phiL1,uL,[],[],[],[],1e-12*ones(1,length(rho1))',[inf*ones(1,length(rho1))]'); %lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
        end

        th(:,ii) = rho1;
        if (ii>1)
            dth = th(:,ii) - th(:,ii-1);
            dth = dth'*dth/length(dth);
        end
        ii=ii+1;
    end
    Ndth = norm(dth);
end
rho = rho1;
end

%% VRFT - Cascade 2
function rho = VRFT_Cascade_MP2(u1,yi1,yi2,ye1,ye2,Cibar,Cebar,Tdi,Tde)

if (isproper(Tdi) ~= 1)
    error('Tdi must be proper');
end

Li = (1-Tdi);
phiLi1 = lsim(Li*Cibar*(1-Tdi),yi1);
phiLi2 = lsim(Li*Cibar*(1-Tdi),yi2);
uLi = lsim(Li*Tdi,u1);
rhoi = (phiLi2'*phiLi1)^-1*(phiLi2'*uLi);
Ciz2 = rhoi'*Cibar;

ri1 = u1 + lsim(Ciz2,yi1);

Le = minreal((1-Tdi)*(1-Tde));
phiLe1 = lsim(Le*Ciz2*Cebar*(1-Tde),ye1);
phiLe2 = lsim(Le*Ciz2*Cebar*(1-Tde),ye2);
uLe = lsim(Le*Tde,ri1);
rhoe = (phiLe2'*phiLe1)^-1*(phiLe2'*uLe);

rho = [rhoe;rhoi];

end
