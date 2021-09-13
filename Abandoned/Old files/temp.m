[ nvar,r1,r2,gammax,b,IB,invAib,cib,free,j ]  = mypreqr3(cov{:,:},y,tau);
[ right2,j ] = multiaddmyqr3(dat{:,2},nvar,r1,r2,gammax,invAib,cib,b,IB,free,y);

[ nvar,r1,r2,gammax,b,IB,invAib,cib,free,j,c,A,right3]  = mypreqr3_alt([dat{:,1} cov{:,:}],y,tau,1);
[ right3,j ] = multialtermyqr3(dat{:,2},nvar,r1,r2,gammax,invAib,A,b,IB,free,y);


%% alter
[ nvar,r1,r2,gammax,b,IB,invAib,cib,free,j,c,A,right3,invAzero]  = mypreqr3_alt([cov{:,:} dat{:,1}],y,tau,1);
j1a=0;j2a=0;
A=sparse(A);invAzero=sparse(invAzero);
tic
for i=2:1000
    [right1,j1,j2,r1,r2,gammax,invAib,A,b,IB,cib ]= multialtermyqr3(dat{:,i},nvar,r1,r2,gammax,invAib,A,b,c,cib,IB,free,y,invAzero);
    j1a=j1a+j1;j2a=j2a+j2;
end
toc

%% add
[ nvar,r1,r2,gammax,b,IB,invAib,cib,free,j ]  = mypreqr3(cov{:,:},y,tau);
tic
ja=0;
for i=1:500
    [ right2,j ] = multiaddmyqr3(dat{:,i},nvar,r1,r2,gammax,invAib,cib,b,IB,free,y);
    ja=ja+j;
end
toc

%% raw
tic
for i=1:nv
    [ right3 ] = myqr3([cov x(:,i)],y,tau);
end
toc

%% inter method
tic
for i=1:nv
    [right4]= rq_fnm([ones(n,1) cov x(:,i)],y,tau);
end
toc

%% simu
%%%%%%%%%%%%
nv = 10000;
rho = 0.75;
n = 200;
sigma = rho*ones(nv,nv);
sigma = sigma+(1-rho)*eye(nv,nv);
mu = 2*ones(1,nv);
rng default  % For reproducibility
x = mvnrnd(mu,sigma,n);

ncov = 20;
cov = normrnd(0,1,n,ncov);

y = normrnd(0,1,n,1);

aibs=sparse(aib);
tic
for i=1:1000
    inv(aibs);
end
toc

