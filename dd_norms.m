function [Norm,varargout] = dd_norms(u,y,Ts,N,pnorm)

data = iddata(y,u,Ts);
% opt = impulseestOptions('pw',5,'RegulKernel','DC');
% sys = impulseest(data,N,opt);

sys = impulseest(data,N);
Mf = sys.Numerator;
% N = length(Mf)
varargout{1} = Mf;

if pnorm == 1
    Norm = sum(abs(Mf));
else if pnorm == 2
    Norm = sqrt(sum(Mf*Mf'));
    else
        S = zeros(N,N);
        for k1 = 1:N
            S(k1:N,k1) = Mf(1:end-k1+1);
        end
        varargout{2} = S;
        Norm = sqrt(max(eig(S'*S)));
    end
end
end