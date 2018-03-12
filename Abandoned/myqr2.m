function [ estimate,j ] = myqr2(x,y,tau)
% beta of X can be negative, improved version



m=length(y);
[~,nvar]=size(x); % all variables in x

% c in min(c^T X)
c=[zeros((1+nvar),1); tau*ones(m,1); (1-tau)*ones(m,1)];

% A matrix & b matrix
A=[ones(m,1) x eye(m,m) -eye(m,m)];
b=y;
A(b<0,:)=-A(b<0,:);
b(b<0)=-b(b<0);

% Ib, Id
IB=[(y>=0).*[(1:m)'+(1+nvar)]+(y<0).*[(1:m)'+(1+nvar)+m]]';
% ID=setdiff(1:(2*(1+nvar)+2*m),IB);

% B matrix
B=A(:,IB); % B=eye(m,m)
cB=c(IB);

TB=[inv(B) zeros(m,1); -cB'/B 1] * [A b; c' 0];
% TB=sparse(TB);

j=1; freevarrow=false(m,1);
while (1)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    rr=r;
    rr(1:1+nvar)=-abs(rr(1:1+nvar));
    if (all(rr>=-1e-7))
        %         sprintf('done')
        break
    end
    % step 3
    % t=find(r==min(r(ID)),1);
    [~,t]=min(rr);
    
    if r(t)<0
        % step 4
        yy=TB(1:(end-1),t);
        if (all(yy<=1e-7))
            sprintf('unboundness')
            break
        end
        % step 5
        b=TB(1:(end-1),end);
        k=b./yy;
        %         k(yy<=1e-7 | ismember(IB,1:nvar+1)' )=inf;
        k=find(k==min(k(yy>0 & ~freevarrow ))       ,1);
        %          [~,k]=min(k);
        %         k=find(k==min(k(IB~=1 & IB~=2)),1);
    else
        % step 4
        yy=TB(1:(end-1),t);
        if (all(yy>=-1e-7))
            sprintf('unboundness')
            break
        end
        % step 5
        b=TB(1:(end-1),end);
        k=-b./yy;
        %         k(yy>=-1e-7 | ismember(IB,1:nvar+1)')=inf;
        k=find(k==min(k(yy<0 & ~freevarrow ))       ,1);
        %          [~,k]=min(k);
        %         k=find(k==min(k(IB~=1 & IB~=2)),1);
    end
    %     [~,k]=min(k);
    % k=find(k==min(k(k>=0)),1);
    
    % step 6
    E=speye(m+1,m+1);
    E(:,k)=[-TB(:,t)./TB(k,t)];
    E(k,k)=1/TB(k,t);
    
    if t<=1+nvar
        freevarrow(k)=true;
    end
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    TB=E * TB;
    %     TB(abs(TB)<1e-8)=0; % too time consuming
    j=j+1;
end

b=TB(1:(end-1),end);
estimate=zeros((nvar+1),1);
for i=1:(nvar+1)
    if isempty(b(IB==i))
        continue
    end
    estimate(i)=b(IB==i);
end
end

