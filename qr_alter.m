function [estimate, pvalue] = qr_alter(z,x,y,tau,nchg,varargin)
% function [ nvar,r1,r2,gammax,b,IB,invAib,cib,freevarrow, j,c,A,estimate, invAzero ] = mypreqr3_alt(x,y,tau,nalt)
% beta of X can be negative, improved version
% sequence: intercept beta beta+ u beta- v

validateattributes(z, {'numeric'},{'2d'},'qr_alter','covariable z',1)
validateattributes(x, {'numeric'},{'2d'},'qr_alter','variable x',2)
[m,nx] = size(x);
validateattributes(y, {'numeric'},{'size',[m,1]},'qr_add','response var y',3)
validateattributes(tau, {'numeric'}, {'scalar','>',0,'<',1},'qr_add','quantile level tau',4)
validateattributes(nchg, {'numeric'}, {'scalar','>',0},'qr_add','number of changed variables',5)
rounds = nx/nchg;
if round(rounds)~=rounds
    error('Number of columns of x is not a multiple of ''nchg''.')
end

inp = inputParser;
addParameter(inp, 'test', [], @(x)any(validatestring(x,{'kernel','wald'})))
addParameter(inp,'hs',true,@islogical)
addParameter(inp,'tol',1e-14,@isscalar)
addParameter(inp,'maxit',10000,@isscalar)
parse(inp, varargin{:})

xback = x;
x = [z x(:,1:nchg)];
[~,nvar]=size(x); % all variables in x

% c in min(c^T X)
c=[zeros(1+nvar,1); tau*ones(m,1); zeros(nchg,1); (1-tau)*ones(m,1)];

% A matrix & b matrix
gammax=[ones(m,1) x];

b=y;
gammax(b<0,:)=-gammax(b<0,:);
b(b<0)=-b(b<0);

% Ib, Id
IB=[(y>=0).*[(1:m)'+1+nvar]+(y<0).*[(1:m)'+1+nvar+m+nchg]];
% ID=setdiff(1:(2*(1+nvar)+2*m),IB);

% TB=[inv(B) zeros(m,1); -cB'/B 1] * [A b; c' 0];

gammax=[gammax; -c(IB')'*gammax];
freevarrow=false(m,1);
% r=[1:nvar+1; zeros(1,nvar+1)];
r1=[1:nvar+1];
r2=[zeros(1,nvar+1-nchg), -((1:nchg)+1+nvar+m)];

j=0;
while (j<inp.Results.maxit)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-inp.Results.tol)
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
            % step 5
            k=b./yy;
            k=find(k==min(k(yy>0 & ~freevarrow )) & (yy>0) ,1);
        else
            % step 5
            k=-b./yy;
            k=find(k==min(k(yy<0 & ~freevarrow )) & (yy<0),1);
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
        % step 5
        k=b./yy;
        k=find(k==min(k(yy>0 & ~freevarrow )) & (yy>0),1);
        
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
        r1(t_rr)=IB(k);
        if IB(k)<=nvar+1-nchg || IB(k)>1+nvar
            r2(t_rr)=IB(k)+m+nchg;
        else
            r2(t_rr)=-(IB(k)+m+nchg);
        end
    else
        if IB(k)<=1+nvar+m+nchg
            gammax(:,t_rr)=full(sparse(k,1,-1,m+1,1));
            r2(t_rr)=-IB(k);
        else
            gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
            r2(t_rr)=IB(k);
        end
        r1(t_rr)=IB(k)-m-nchg;
    end
    
    gammax=gammax-ee*gammax(k,:);
    b=b-ee(1:end-1)*b(k);
    IB(k)=t;
    j=j+1;
end

if j==inp.Results.maxit
    warning('May not converge.')
end

cib=c(IB);
A=[ones(m,1) x eye(m,m) -x(:,end-nchg+1:end) -eye(m,m)];
A(y<0,:)=-A(y<0,:);
invAib = inv(A(:,IB));
estimate=zeros((nvar+1),1);

for i=1:(nvar+1-nchg)
    try
        estimate(i)=b(IB==i);
    catch
        estimate(i)=0;
        sprintf('Warning!')
    end
end
for i=nvar+2-nchg:nvar+1
    try
        if isempty(b(IB==i))
            estimate(i)=-b(IB==(i+m+nchg));
        else
            estimate(i)=b(IB==i);
        end
    catch
        estimate(i)=0;
        sprintf('Warning!')
    end
    
end
invAzero = [zeros(1+nvar,m-1-nvar); eye(m-1-nvar)];

bd = band(tau, m, inp.Results.hs);
h1 = (norminv(tau + bd) - norminv(tau - bd));

pvalue = [];
for i = nchg+1:nchg:(nx-nchg+1)
%     i
    estimate = [estimate multialtermyqr3(xback(:,i:i+nchg-1),nvar,r1,r2,...
        gammax,invAib,A,b,c,cib,IB,freevarrow,y,invAzero,inp.Results.tol,inp.Results.maxit)];
end
if ~isempty(inp.Results.test)
    esti = 1;
    switch inp.Results.test
        case 'kernel'
            for i = 1:nchg:(nx-nchg+1)
                pvalue = [pvalue ker_add(z, xback(:,i:(i+nchg-1)), y, estimate(:,esti), tau, h1)'];
                esti = esti+1;
            end
        case 'wald'
            for i = 1:nchg:(nx-nchg+1)
                pvalue = [pvalue wald_add(z, xback(:,i:(i+nchg-1)), y, estimate(:,esti), bd, tau)'];
                esti = esti+1;
            end
    end
end
end

function [ estimate,j1,j2,r1,r2,gammax,invAib,A,b,IB,cib ] = multialtermyqr3(xnew,nvar,r1,r2,gammax,...
    invAib,A,b,c,cib,IB,freevarrow,y,invAzero,tol,maxit)
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
while (j1<maxit)
    % start iteration
    % step 2
    rr=gammax(end,:);
    %     rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    rr=[rr; -rr];
    %     rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-tol)
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
    k=find(k==min(k(yy>0 & ~freevarrow )) & (yy>0),1);
    
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


if j1==maxit
    warning('May not converge.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j2=0;
while (j2<maxit)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2>0)-rr.*(r2<0)];
    %     rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-tol)
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
    k=find(k==min(k(yy>0 & ~freevarrow )) & (yy>0),1);
    
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


if j2==maxit
    warning('May not converge.')
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
invAib=invA(A(:,IB),IB,m,y,nvar,nalt,invAzero);
% invAib = inv(A(:,IB));
cib = c(IB);

end

function invA = invA(Aib,IB,m,y,nvar,nalt,invAzero)
% x adds ones,  if y<0 x=-x
ys=y<0;

R=find(IB <= (1+nvar) | (IB > 1+nvar+m & IB <= 1+nvar+m+nalt));
IB(R)=[];

IBs = IB> 1+nvar+m;
R=[R' my_setdiff(1:m,R)];

IB(IBs)=IB(IBs)-m-nalt;
IB=IB-1-nvar;
ysIB=ys(IB);

negsign =[false(1+nvar,1); xor(IBs,ysIB)];
IB = [my_setdiff(1:m,IB) IB'];

x12=Aib(IB(1:1+nvar),R(1:1+nvar));
x21=Aib(IB(2+nvar:end),R(1:1+nvar));

b112=inv(x12);
invAtmp=[[b112; -x21*b112] invAzero];
invAtmp(negsign,:)=-invAtmp(negsign,:);
[~,IBsinv]=sort(IB);
[~,Rsinv]=sort(R);
invA=invAtmp(Rsinv,IBsinv);
end

function C = my_setdiff(A,B)
    bits = false(1, max(max(A), max(B)));
    bits(A) = true;
    bits(B) = false;
    C = A(bits(A));
end

function [ pvalue ] = ker_add(x, xnew, y, est, tau, h1)
p = length(est);
n = size(x,1);
rdf = n - p;

xx = [ones(n,1) x xnew];
uhat = y - xx * est;

h = h1 .* min(sqrt(var(uhat)), (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34);
f = uhat./h;
f = exp(-0.5 * f.^2) ;
f = f./ sqrt(2*pi);
f = f./h;

fxxinv = bsxfun(@times,sqrt(f),xx);
fxxt = fxxinv';
fxxinv = inv(fxxt * fxxinv);

cov = tau .* (1 - tau) .* fxxinv * xx' * xx * fxxinv;

serr = sqrt(diag(cov));
tvalue = -abs(est./serr);

p = NaN(size(tvalue)); % single if p or v is
xsq = tvalue.^2;
t = (rdf < xsq);
if any(t(:))
    pt(t) = betainc(rdf ./ (rdf + xsq(t)), rdf/2, 0.5, 'lower') / 2;
end
t = ~t;
if any(t(:))
    pt(t) = betainc(xsq(t) ./ (rdf + xsq(t)), 0.5, rdf/2, 'upper') / 2;
end
pvalue = 2 * pt;
end

function [ pvalue ] = wald_add(x, xnew, y, est, bandwidth, tau)
p = length(est);
n = size(x,1);
xx = [ones(n,1) x xnew];
xxt = xx';
xxi = xxt * xx;
resid = y - xx * est;
rdf = n - p;

pz = sum(abs(resid) < sqrt(eps));
h = max(p + 1, ceil(n * bandwidth));
ir = (pz + 1):(h + pz + 1);
[~,Id]=sort(abs(resid));
Id = resid(Id);
ord_resid = sort(Id(ir));

xt = ir'./(n - p);
% sparsity = myqr3(xt', ord_resid, tau);
sparsity = rq_fnm([ones(h+1,1) xt], ord_resid, 0.5);
sparsity = sparsity(2);
cov = sparsity.^2 .* tau .* (1 - tau).* inv(xxi) ;
serr = sqrt(diag(cov));
tvalue = -abs(est./serr);

p = NaN(size(tvalue)); % single if p or v is
xsq = tvalue.^2;
t = (rdf < xsq);
if any(t(:))
    pt(t) = betainc(rdf ./ (rdf + xsq(t)), rdf/2, 0.5, 'lower') / 2;
end
t = ~t;
if any(t(:))
    pt(t) = betainc(xsq(t) ./ (rdf + xsq(t)), 0.5, rdf/2, 'upper') / 2;
end

pvalue = 2 * pt;
end

function b = band(p, m, hs)
x0 = norminv(p);
f0 = normpdf(x0);
if hs
    b=m^(-1/3) * norminv(1 - 0.05/2)^(2/3) * ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3);
else
    b=m^(-0.2) * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^0.2;
end
end

function [b,it] = rq_fnm(X, y, p)
% Construct the dual problem of quantile regression
% Solve it with lp_fnm
%
%
[m n] = size(X);
u = ones(m, 1);
a = (1 - p) .* u;
[b,it] = lp_fnm(X', -y', X' * a, u, a);
b = -b';
end

function [y,it] = lp_fnm(A, c, b, u, x)
% Solve a linear program by the interior point method:
% min(c * u), s.t. A * x = b and 0 < x < u
% An initial feasible solution has to be provided as x
%
% Function lp_fnm of Daniel Morillo & Roger Koenker
% Translated from Ox to Matlab by Paul Eilers 1999
% Modified by Roger Koenker 2000--
% More changes by Paul Eilers 2004


% Set some constants
beta = 0.9995;
small = 1e-8;
max_it = 500;
[m n] = size(A);

% Generate inital feasible point
s = u - x;
y = (A' \  c')';
r = c - y * A;
r = r + 0.001 * (r == 0);    % PE 2004
z = r .* (r > 0);
w = z - r;
gap = c * x - y * b + w * u;

% Start iterations
it = 0;
while (gap) > small & it < max_it
    it = it + 1;
    
    %   Compute affine step
    q = 1 ./ (z' ./ x + w' ./ s);
    r = z - w;
    Q = spdiags(sqrt(q), 0, n, n);
    AQ = A * Q;          % PE 2004
    rhs = Q * r';        % "
    dy = (AQ' \ rhs)';   % "
    dx = q .* (dy * A - r)';
    ds = -dx;
    dz = -z .* (1 + dx ./ x)';
    dw = -w .* (1 + ds ./ s)';
    
    % Compute maximum allowable step lengths
    fx = bound(x, dx);
    fs = bound(s, ds);
    fw = bound(w, dw);
    fz = bound(z, dz);
    fp = min(fx, fs);
    fd = min(fw, fz);
    fp = min(min(beta * fp), 1);
    fd = min(min(beta * fd), 1);
    
    % If full step is feasible, take it. Otherwise modify it
    if min(fp, fd) < 1
        
        % Update mu
        mu = z * x + w * s;
        g = (z + fd * dz) * (x + fp * dx) + (w + fd * dw) * (s + fp * ds);
        mu = mu * (g / mu) ^3 / ( 2* n);
        
        % Compute modified step
        dxdz = dx .* dz';
        dsdw = ds .* dw';
        xinv = 1 ./ x;
        sinv = 1 ./ s;
        xi = mu * (xinv - sinv);
        rhs = rhs + Q * (dxdz - dsdw - xi);
        dy = (AQ' \ rhs)';
        dx = q .* (A' * dy' + xi - r' -dxdz + dsdw);
        ds = -dx;
        dz = mu * xinv' - z - xinv' .* z .* dx' - dxdz';
        dw = mu * sinv' - w - sinv' .* w .* ds' - dsdw';
        
        % Compute maximum allowable step lengths
        fx = bound(x, dx);
        fs = bound(s, ds);
        fw = bound(w, dw);
        fz = bound(z, dz);
        fp = min(fx, fs);
        fd = min(fw, fz);
        fp = min(min(beta * fp), 1);
        fd = min(min(beta * fd), 1);
        
    end
    
    % Take the step
    x = x + fp * dx;
    s = s + fp * ds;
    y = y + fd * dy;
    w = w + fd * dw;
    z = z + fd * dz;
    gap = c * x - y * b + w * u;
    %disp(gap);
end
end

function b = bound(x, dx)
% Fill vector with allowed step lengths
% Support function for lp_fnm
b = 1e20 + 0 * x;
f = find(dx < 0);
b(f) = -x(f) ./ dx(f);
end