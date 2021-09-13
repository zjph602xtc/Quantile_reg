function [estimate, pvalue] = qr_add(z,x,y,tau,nchg,varargin)
% beta of X can be negative, improved version

validateattributes(z, {'numeric'},{'2d'},'qr_add','covariable z',1)
validateattributes(x, {'numeric'},{'2d'},'qr_add','variable x',2)
[~,nx] = size(x);
[m,nvar]=size(z);
validateattributes(y, {'numeric'},{'size',[m,1]},'qr_add','response var y',3)
validateattributes(tau, {'numeric'}, {'scalar','>',0,'<',1},'qr_add','quantile level tau',4)
validateattributes(nchg, {'numeric'}, {'scalar','>',0},'qr_add','number of added variables',5)
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

% c in min(c^T X)
c=[zeros((1+nvar),1); tau*ones(m,1); (1-tau)*ones(m,1)];

% A matrix & b matrix
gammax=[ones(m,1) z];
b=y;
gammax(b<0,:)=-gammax(b<0,:);
b(b<0)=-b(b<0);

% Ib, Id
IB=[(y>=0).*[(1:m)'+(1+nvar)]+(y<0).*[(1:m)'+(1+nvar)+m]];
% ID=setdiff(1:(2*(1+nvar)+2*m),IB);

% TB=[inv(B) zeros(m,1); -cB'/B 1] * [A b; c' 0];

gammax=[gammax; -c(IB')'*gammax];
freevarrow=false(m,1);
% r=[1:nvar+1; zeros(1,nvar+1)];
r1=1:nvar+1;
r2=zeros(1,nvar+1);

j=0;
while (j<inp.Results.maxit)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2~=0)];
    rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-inp.Results.tol)
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
        ga=gammax(1:end-1,t_rr);
        if gammax(end,t_rr)<0
            % step 5
            k=b./ga;
            k=find(k==min(k(ga>0 & ~freevarrow )) & (ga>0) ,1);
        else
            % step 5
            k=-b./ga;
            k=find(k==min(k(ga<0 & ~freevarrow )) & (ga<0) ,1);
        end
        % pivoting
        % step 6
        ga=gammax(:,t_rr);
        freevarrow(k)=true;
    else % r(2,t_rr)~=0
        % choose k
        % step 4
        if tsep==1
            ga=gammax(1:end-1,t_rr);
        else
            ga=-gammax(1:end-1,t_rr);
        end
        % step 5
        k=b./ga;
        k=find(k==min(k(ga>0 & ~freevarrow )) & (ga>0),1);
        
        % pivoting
        % step 6
        if tsep==1
            ga=gammax(:,t_rr);
        else
            ga=[ga; 1-gammax(end,t_rr)];
        end
        
    end
    
    ee=ga./ga(k);
    ee(k)=1-1/ga(k); % that is - ee.
    
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

if j==inp.Results.maxit
    warning('May not converge.')
end

cib=c(IB);
invAib=[ones(m,1) z eye(m,m) -eye(m,m)];
invAib(y<0,:)=-invAib(y<0,:);
invAib=invAib(:,IB);
% invAib = inv(invAib);

bd = band(tau, m, inp.Results.hs);
h1 = (norminv(tau + bd) - norminv(tau - bd));

estimate = [];
pvalue = [];
x(y<0,:)=-x(y<0,:);
% BinvX=invAib\x; % slower
invAib = inv(invAib);
BinvX = invAib*x;
for i = 1:nchg:(nx-nchg+1)
    estimate = [estimate multiaddmyqr3(x(:,i:(i+nchg-1)),nvar,r1,r2,gammax,...
        BinvX(:,i:(i+nchg-1)),cib,b,IB,freevarrow,y,inp.Results.tol,inp.Results.maxit)];
end
if ~isempty(inp.Results.test)
    esti = 1;
    switch inp.Results.test
        case 'kernel'
            for i = 1:nchg:(nx-nchg+1)
                pvalue = [pvalue ker_add(z, x(:,i:(i+nchg-1)), y, estimate(:,esti), tau, h1)'];
                esti = esti+1;
            end
        case 'wald'
            for i = 1:nchg:(nx-nchg+1)
                pvalue = [pvalue wald_add(z, x(:,i:(i+nchg-1)), y, estimate(:,esti), bd, tau)'];
                esti = esti+1;
            end
    end
end

end

function [ estimate,j ] = multiaddmyqr3(xnew,nvar,r1,r2,gammax,...
    BinvX,cib,b,IB,freevarrow,y,tol,maxit)
% intercept old x  new x
[~,newnvar]=size(xnew);
nvar =nvar+newnvar;
m=length(y);
% xnew(y<0,:)=-xnew(y<0,:);

% BinvX=invAib*xnew;

gammax = [gammax [BinvX; -cib'*BinvX]];
r1 = [r1 2*m+nvar+2:2*m+nvar+1+newnvar];
r2 = [r2 zeros(1, newnvar)];

j=0;
while (j<maxit)
    % start iteration
    % step 2
    rr=gammax(end,:);
    rr=[rr; (1-rr).*(r2~=0)];
    rr(1,r2==0)=-abs(rr(1,r2==0));
    
    if all(rr(:)>=-tol)
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
            % step 5
            k=b./yy;
            k=find(k==min(k(yy>0 & ~freevarrow )) & (yy>0) ,1);
        else
            % step 5
            k=-b./yy;
            k=find(k==min(k(yy<0 & ~freevarrow )) & (yy<0) ,1);
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

if j==maxit
    warning('May not converge.')
end

estimate=[];
for i=[1:nvar-newnvar+1 2*m+nvar+2:2*m+nvar+1+newnvar]
    estimatenew = b(IB==i);
    if isempty(estimatenew)
        estimatenew = 0;
    end
    estimate=[estimate; estimatenew];
end
end

function [ pvalue ] = ker_add(x, xnew, y, est, tau, h1)
p = length(est);
n = size(x,1);
rdf = n - p;

xx = [ones(n,1) x xnew];
uhat = y - xx * est;

h = h1 .* min(sqrt(var(uhat)), iqr(uhat)/1.34);
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