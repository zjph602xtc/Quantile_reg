function [res] = mycompositerq(dat,tau,testvar,response,covariate)
%function 
% testvar='x';
% response='y';
% covariate={'z'};
%%%%%
% response first testvar second  covariate finally
n=height(dat);
nvar=2+length(covariate);
ntau=length(tau);
objectin=zeros(1,(2*n+nvar-1)*ntau+1);
betaind=[];
for i=1:ntau
    objectin(((i-1)*(2*n+nvar-1)+1+nvar):(n+(i-1)*(2*n+nvar-1)+nvar)) = tau(i);
    objectin((n+(i-1)*(2*n+nvar-1)+nvar+1):(2*n+(i-1)*(2*n+nvar-1)+nvar)) = 1-tau(i);
    betaind = [betaind, (2+(i-1)*(nvar-1+2*n)):((i-1)*(nvar-1+2*n)+nvar)];
end
betaind = [1,betaind];
%%%%%constraints
if isempty(covariate)
    designx=ones(n,1);
else
    designx=[ones(n,1), dat{:,covariate}];
end
Xd=dat{:,testvar};
Id=speye(n);
dtem=[designx,Id,-Id];
dtemp=repmat({dtem},1,ntau);
constmat=[repmat(Xd,ntau,1) sparse(blkdiag(dtemp{:}))];

constrhs=repmat(dat{:,response},ntau,1);
I=[1:(2*n+nvar-1)*ntau+1];
I(betaind)=[];
lower=repmat(-inf,1,(2*n+nvar-1)*ntau+1);
lower(I)=0;

res=linprog(objectin,[],[],constmat,constrhs,lower,[],[],optimoptions('linprog','Algorithm','dual-simplex','Display','off'));
res=res(betaind);
end

