function [ estimate,j ] = multiaddmyqr3(xnew,nvar,r1,r2,gammax,invAib,cib,b,IB,freevarrow,y)

[~,newnvar]=size(xnew);
nvar =nvar+newnvar;
m=length(y);
xnew(y<0,:)=-xnew(y<0,:);

BinvX=invAib*xnew;


gammax = [gammax [BinvX; -cib'*BinvX]];
r1 = [r1 2*m+nvar+2:2*m+nvar+1+newnvar];
r2 = [r2 zeros(1, newnvar)];



j=0;
while (1)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2~=0)];
    rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-1e-7)
        %         sprintf('done')
        break
    end
    % step 3
    [tsep,t_rr]=find(rr==min(rr(:)),1);
    %     t=r(tsep,t_rr);
    if tsep ==1
        t=r1(t_rr);
    else
        t=r2(t_rr);
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
        else
            yy=[yy; 1-gammax(end,t_rr)];
        end
        
    end
    
    ee=yy./yy(k);
    ee(k)=1-1/yy(k); % that is - ee.
    
    if IB(k)<=nvar+m+1
        gammax(:,t_rr)=full(sparse(k,1,1,m+1,1));
        %         r(:,t_rr)=[IB(k); IB(k)+m];
        r1(t_rr)=IB(k);
        r2(t_rr)=IB(k)+m;
    else
        gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
        %         r(:,t_rr)=[IB(k)-m; IB(k)];
        r1(t_rr)=IB(k)-m;
        r2(t_rr)=IB(k);
    end
    
    gammax=gammax-ee*gammax(k,:);
    b=b-ee(1:end-1)*b(k);
    IB(k)=t;
    
    j=j+1;
end

estimate=[];
for i=[1:nvar 2*m+nvar+2:2*m+nvar+1+newnvar]
    %     if isempty(b(IB==i))
    %         continue
    %     end
    estimate=[estimate; b(IB==i)];
end
end

