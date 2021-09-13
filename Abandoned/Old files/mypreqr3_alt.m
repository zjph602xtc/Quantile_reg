function [ nvar,r1,r2,gammax,b,IB,invAib,cib,freevarrow, j,c,A,estimate, invAzero ] = mypreqr3_alt(x,y,tau,nalt)
% beta of X can be negative, improved version
% sequence: intercept beta beta+ u beta- v


m=length(y);
[~,nvar]=size(x); % all variables in x

% c in min(c^T X)
c=[zeros(1+nvar,1); tau*ones(m,1); zeros(nalt,1); (1-tau)*ones(m,1)];

% A matrix & b matrix
gammax=[ones(m,1) x];

b=y;
gammax(b<0,:)=-gammax(b<0,:);
b(b<0)=-b(b<0);

% Ib, Id
IB=[(y>=0).*[(1:m)'+1+nvar]+(y<0).*[(1:m)'+1+nvar+m+nalt]];
% ID=setdiff(1:(2*(1+nvar)+2*m),IB);

% TB=[inv(B) zeros(m,1); -cB'/B 1] * [A b; c' 0];

gammax=[gammax; -c(IB')'*gammax];
freevarrow=false(m,1);
% r=[1:nvar+1; zeros(1,nvar+1)];
r1=[1:nvar+1];
r2=[zeros(1,nvar+1-nalt) -((1:nalt)+1+nvar+m)];

j=1;
while (1)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-1e-7)
        %         sprintf('done')
        break
    end
    % step 3
    [tsep,t_rr]=find(rr==min(rr(:)),1);
    %     t=r(tsep,t_rr);
    if tsep==1
        t=r1(t_rr);
    else
        t=abs(r2(t_rr));
    end
    
    
    if r2(t_rr)==0
        % choose k
        % step 4
        yy=gammax(1:end-1,t_rr);
        if gammax(end,t_rr)<0
            %             if (all(yy<=1e-7))
            %                 sprintf('unboundness')
            %                 break
            %             end
            % step 5
            k=b./yy;
            k=find(k==min(k(yy>0 & ~freevarrow )) ,1);
        else
            % step 4
            %             if (all(yy>=-1e-7))
            %                 sprintf('unboundness')
            %                 break
            %             end
            % step 5
            k=-b./yy;
            k=find(k==min(k(yy<0 & ~freevarrow )) ,1);
        end
        
        % pivoting
        % step 6
        yy=gammax(:,t_rr);
        freevarrow(k)=true;
        
    else % r(2,t_rr)~=0
        % choose k
        % step 4
        if tsep==1
            yy=gammax(1:end-1,t_rr);
        else
            yy=-gammax(1:end-1,t_rr);
        end
        %         if (all(yy<=1e-7))
        %             sprintf('unboundness')
        %             break
        %         end
        % step 5
        k=b./yy;
        k=find(k==min(k(yy>0 & ~freevarrow )),1);
        
        % pivoting
        % step 6
        if tsep==1
            yy=gammax(:,t_rr);
        elseif r2(t_rr)>0
            yy=[yy; 1-gammax(end,t_rr)];
        else
            yy=[yy; -gammax(end,t_rr)];
        end
        
    end
    
    ee=yy./yy(k);
    ee(k)=1-1/yy(k); % that is - ee.
    
    if IB(k)<=nvar+m+1
        gammax(:,t_rr)=full(sparse(k,1,1,m+1,1));
        %         r(:,t_rr)=[IB(k); IB(k)+m];
        r1(t_rr)=IB(k);
        if IB(k)<=nvar+1-nalt || IB(k)>1+nvar
            r2(t_rr)=IB(k)+m+nalt;
        else
            r2(t_rr)=-(IB(k)+m+nalt);
        end
    else
        %         gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
        %         r(:,t_rr)=[IB(k)-m; IB(k)];
        if IB(k)<=1+nvar+m+nalt
            gammax(:,t_rr)=full(sparse(k,1,-1,m+1,1));
            r2(t_rr)=-IB(k);
        else
            gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
            r2(t_rr)=IB(k);
        end
        r1(t_rr)=IB(k)-m-nalt;
    end
    
    gammax=gammax-ee*gammax(k,:);
    b=b-ee(1:end-1)*b(k);
    IB(k)=t;
    j=j+1;
end
cib=c(IB);
A=[ones(m,1) x eye(m,m) -x(:,end-nalt+1:end) -eye(m,m)];
A(y<0,:)=-A(y<0,:);
invAib = inv(A(:,IB));


estimate=zeros((nvar+1),1);

for i=1:(nvar+1-nalt)
    %     if isempty(b(IB==i))
    %         continue
    %     end
    try
        estimate(i)=b(IB==i);
    catch
        estimate(i)=0;
        sprintf('Warning!')
    end
end
for i=nvar+2-nalt:nvar+1
    try
        if isempty(b(IB==i))
            estimate(i)=-b(IB==(i+m+nalt));
        else
            estimate(i)=b(IB==i);
        end
    catch
        estimate(i)=0;
        sprintf('Warning!')
    end
    
end
invAzero = [zeros(1+nvar,m-1-nvar); eye(m-1-nvar)];
end

