function [ estimate,j1,j2,r1,r2,gammax,invAib,A,b,IB,cib ] = multialtermyqr3(xnew,nvar,r1,r2,gammax,...
    invAib,A,b,c,cib,IB,freevarrow,y,invAzero)
% change last nalt columns!

[~,nalt]=size(xnew);
m=length(y);
xnew(y<0,:)=-xnew(y<0,:);

r1 = [r1 (1:nalt)+nvar+1-nalt];
r2= [r2 -((1:nalt)+nvar+1+m)];

IBold = intersect(IB, [nvar-nalt+2:nvar+1 (nvar-nalt+2:nvar+1)+m+nalt]); % this is not real IB
% IDold = setdiff([nvar+1:nvar+nalt (nvar+1:nvar+nalt)+m+nalt], IB);
naltpool=(1:nalt)+1+nvar+2*m+nalt;
for jj = 1:length(naltpool)
    IB(IB==IBold(jj)) = naltpool(jj);
end

A(:,[nvar+2-nalt:nvar+1 (nvar+2-nalt:nvar+1)+m+nalt])=[xnew -xnew];
gammax = [gammax [invAib*xnew; zeros(1,nalt)]];

% rr = c(r1,:)'-cib'*invAib*A(:,r1);
cstar = [zeros(1,1+nvar+nalt+2*m) ones(1,nalt)];
rstar = -cstar(IB)*invAib*A(:,r1);
gammax = [gammax;  rstar];

rrr = [c'-cib'*invAib*A];
gammax(end-1,:)=rrr(:,r1);

j1=0; naltdet = 0;
while (1)
    % start iteration
    % step 2
    rr=gammax(end,:);
    %     rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    rr=[rr; -rr];
    %     rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-1e-7)
        %         sprintf('done')
        break
    end
    % step 3
    [tsep,t_rr]=find(rr==min(rr(:)),1);
    if tsep==1
        t=r1(t_rr);
    else
        t=abs(r2(t_rr));
    end
    
    
    if tsep==1
        yy=gammax(1:end-2,t_rr);
    else
        yy=-gammax(1:end-2,t_rr);
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
        yy=[yy; 1-gammax(end-1,t_rr)];
        yy=[yy; -gammax(end,t_rr)];
    else
        yy=[yy; -gammax(end-1:end,t_rr)];
    end
    
    %     end
    
    ee=yy./yy(k);
    ee(k)=1-1/yy(k); % that is - ee.
    
    if IB(k)<=nvar+m+1
        gammax(:,t_rr)=full(sparse(k,1,1,m+2,1));
        %         r(:,t_rr)=[IB(k); IB(k)+m];
        r1(t_rr)=IB(k);
        if IB(k)<=nvar+1-nalt || IB(k)>1+nvar
            r2(t_rr)=IB(k)+m+nalt;
        else
            r2(t_rr)=-(IB(k)+m+nalt);
        end
    elseif IB(k)<=nvar+2*m+nalt+1
        if IB(k)<=1+nvar+m+nalt
            gammax(:,t_rr)=full(sparse(k,1,-1,m+2,1));
            r2(t_rr)=-IB(k);
        else
            gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+2,1));
            r2(t_rr)=IB(k);
        end
        r1(t_rr)=IB(k)-m-nalt;
    else
        naltdet = naltdet +1;
        gammax(:,t_rr)=[];
        r1(t_rr)=[];
        r2(t_rr)=[];
        if naltdet==nalt
            gammax(end,:)=[];
            gammax=gammax-ee(1:end-1)*gammax(k,:);
            b=b-ee(1:end-2)*b(k);
            IB(k)=t;
            j1=j1+1;
            break;
        end
    end
    gammax=gammax-ee*gammax(k,:);
    b=b-ee(1:end-2)*b(k);
    IB(k)=t;
    
    j1=j1+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j2=0;
while (1)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    %     rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-1e-7)
        %         sprintf('done')
        break
    end
    % step 3
    [tsep,t_rr]=find(rr==min(rr(:)),1);
    if tsep ==1
        t=r1(t_rr);
    else
        t=abs(r2(t_rr));
    end
    
    % choose k
    % step 4
    if tsep==1
        yy=gammax(1:end-1,t_rr);
    else
        yy=-gammax(1:end-1,t_rr);
    end
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
        gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
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
    
    j2=j2+1;
end


estimate=zeros((nvar+1),1);
for i=1:(nvar+1-nalt)
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
% invAib = inv(A(:,IB));
invAib=invA(A(:,IB),IB,m,y,nvar,nalt,invAzero);
cib = c(IB);
end


