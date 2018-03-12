function [ TB,IB,estimate,j ] = paramyqr(x,y,tau)
% calculate for different tau
m=length(y);
[~,nvar]=size(x); % all variables in x

% c in min(c^T X)
c=[zeros(2*(1+nvar),1); 0*ones(m,1); 1*ones(m,1)];

% A matrix & b matrix
DDx=zeros(m,2*(1+nvar));
DDx(:,1:2:end)=[ones(m,1) x];
DDx(:,2:2:end)=-[ones(m,1) x];
A=[DDx eye(m,m) -eye(m,m)];
b=y;
A(b<0,:)=-A(b<0,:);
b(b<0)=-b(b<0);

% Ib, Id
IB=[find(y>=0)+2*(1+nvar); find(y<0)+2*(1+nvar)+m]';
% ID=setdiff(1:(2*(1+nvar)+2*m),IB);

% B matrix
B=A(:,IB); % B=eye(m,m)
cB=c(IB);

TB=[inv(B) zeros(m,1); -cB'/B 1] * [A b; c' 0];
TB=sparse(TB);

j=1;
while (1)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    if (all(r>=-1e-7))
        sprintf('done')
        break
    end
    % step 3
    % t=find(r==min(r(ID)),1);
    [~,t]=min(r);
    % step 4
    yy=TB(1:(end-1),t);
    if (all(yy<=0))
        sprintf('unboundness')
        break
    end
    % step 5
    b=TB(1:(end-1),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    % k=find(k==min(k(k>=0)),1);
    
    % step 6
    E=speye(m+1,m+1);
    E(:,k)=[-TB(:,t)./TB(k,t)];
    E(k,k)=1/TB(k,t);
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    TB=E * TB;
    % TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
end

b=TB(1:(end-1),end);
estimate=zeros((nvar+1),1);
for i=1:2*(nvar+1)
    if isempty(b(IB==i))
        continue
    end
    if (mod(i,2)==1)
        estimate((i+1)/2)=b(IB==i);
    else
        estimate(i/2)=-b(IB==i);
    end
end

end

