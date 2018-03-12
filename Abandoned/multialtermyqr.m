function [ estimate,j, A,TB,IB ] = multialtermyqr(xnew,A,TB,IB,c,nvar,y)

[~,ntb]=size(TB);
[~,nnew] = size(xnew);
m=length(y);
xnew(y<0,:)=-xnew(y<0,:);
bi=inv(A(:,IB));

% aix = zeros(m, nnew*2);
% aix(:,1:2:end) = xnew;
% aix(:,2:2:end) = -xnew;


bix=bi*xnew;
cib=c(IB);

nobasic=setdiff(3:2+2*nnew,IB);
nobasicodd=nobasic(mod(nobasic,2)~=0);
nobasiceven=nobasic(mod(nobasic,2)==0);

TB(:,nobasicodd)=[bix(:,(nobasicodd-1)/2); zeros(1,length(nobasicodd)) ];
TB(:,nobasiceven)=[-bix(:,(nobasiceven-2)/2); zeros(1,length(nobasiceven)) ];

basic = setdiff(3:2+2*nnew,nobasic);
TB=TB(:,[1:end-1 basic end]);
for ibbb=1:length(basic)
    IB(IB==basic(ibbb))=ntb-1+ibbb;
end
basicodd=basic(mod(basic,2)~=0);
basiceven=basic(mod(basic,2)==0);
TB(:,basicodd)=[bix(:,(basicodd-1)/2); zeros(1,length(basicodd)) ];
TB(:,basiceven)=[-bix(:,(basiceven-2)/2); zeros(1,length(basiceven)) ];

ttt=zeros(m,nnew*2);
ttt(:,1:2:end)=xnew;
ttt(:,2:2:end)=-xnew;
A(:,3:2+2*nnew)=ttt;

cc=[zeros(ntb-1,1); ones(length(basic),1)];
TB(end,[1:ntb-1])=-cc(IB)'*bi*A;
%!

% c=[c ;zeros(length(basic),1)];
% cib=c(IB);

TB=[TB; c'-cib'*bi*A zeros(1,length(basic)+1)];



j=1;
while (j<5000)
    % start iteration
    % step 2
    r=TB(end-1,1:(end-1));
    if (all(r>=-1e-7))
        %         sprintf('done')
        break
    end
    % step 3
    [~,t]=min(r);
%     t=find(r<0,1);
    % step 4
    yy=TB(1:(end-2),t);
    if (all(yy<=0))
        sprintf('unboundness')
        break
    end
    % step 5
    b=TB(1:(end-2),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    %     TB=E * TB;
    TB=TB+ee*TB(k,:);
    
    
    %     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
    if (j==4800); sprintf('not converge stage 1'); end;
end

TB=TB([1:(end-2) end],[1:ntb-1 end]);

while (j<5000)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    if (all(r>=-1e-7))
        %         sprintf('done')
        break
    end
    % step 3
    [~,t]=min(r);
    % step 4
    yy=TB(1:(end-1),t);
    %     if (all(yy<=0))
    %         sprintf('unboundness')
    %         break
    %     end
    % step 5
    b=TB(1:(end-1),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    %     TB=E * TB;
    TB=TB+ee*TB(k,:);
    
    
    %     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
        if (j==4800); sprintf('not converge stage 2'); end;
end


clear estimate;

b=TB(1:(end-1),end);
estimate=[];
estimate1=[];
betaindex=1:2:2*(nvar+1);
for i = 1:nvar+1
    if any(IB==betaindex(i))
        estimate = [estimate; b(IB==betaindex(i))];
        estimate1 = [estimate1; 0];
    end
    if any(IB==betaindex(i)+1)
        estimate = [estimate; 0];
        estimate1 = [estimate1; b(IB==betaindex(i)+1)];
    end
end

estimate=estimate-estimate1;
end


