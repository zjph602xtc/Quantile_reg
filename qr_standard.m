function [ estimate,pvalue,j ] = qr_standard(x,y,tau,varargin)
% beta of X can be negative, improved version
% intercept first in the response

validateattributes(x, {'numeric'},{'2d'},'qr_standard','variable x',1)
[m,nvar]=size(x);
validateattributes(y, {'numeric'},{'size',[m,1]},'qr_standard','response var y',2)
validateattributes(tau, {'numeric'}, {'scalar','>',0,'<',1},'qr_standard','quantile level tau',3)

inp = inputParser;
addParameter(inp, 'test', [], @(x)any(validatestring(x,{'kernel','wald'})))
addParameter(inp,'hs',true,@islogical)
addParameter(inp,'tol',1e-14,@isscalar)
addParameter(inp,'maxit',10000,@isscalar)
addParameter(inp,'method',[],@(x)any(validatestring(x,{'interior'})))

parse(inp, varargin{:})

if isempty(inp.Results.method)
    % c in min(c^T X)
    c=[zeros((1+nvar),1); tau*ones(m,1); (1-tau)*ones(m,1)];
    
    % A matrix & b matrix
    gammax=[ones(m,1) x];
    b=y;
    gammax(b<0,:)=-gammax(b<0,:);
    b(b<0)=-b(b<0);
    
    % Ib, Id
    IB=[(y>=0).*[(1:m)'+(1+nvar)]+(y<0).*[(1:m)'+(1+nvar)+m]];
    
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
            else
                yy=[yy; 1-gammax(end,t_rr)];
            end
            
        end
        
        ee=yy./yy(k);
        ee(k)=1-1/yy(k); % that is - ee.
        
        if IB(k)<=nvar+m+1
            gammax(:,t_rr)=full(sparse(k,1,1,m+1,1));
            r1(t_rr)=IB(k);
            r2(t_rr)=IB(k)+m;
        else
            gammax(:,t_rr)=full(sparse([k m+1],1,[-1 1],m+1,1));
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
    
    estimate=zeros((nvar+1),1);
    for i=1:(nvar+1)
        %     if isempty(b(IB==i))
        %         continue
        %     end
        try
            estimate(i)=b(IB==i);
        catch
            estimate(i)=0;
        end
    end
else
    estimate = rq_fnm([ones(m,1) x],y,tau);
end

% test part
p = 1+nvar;
xx = [ones(m,1) x];
rdf = m - p;
resid = y - xx * estimate;
if ~isempty(inp.Results.test)
    switch inp.Results.test
        case 'kernel'
            h = band(tau, m, inp.Results.hs);
            if (tau+h)>1 || (tau-h)<0
                error('Data do not meet the assumption of kernel test')
            end
            h = (norminv(tau+h)-norminv(tau-h)) .* ...
                min(sqrt(var(resid)), (quantile(resid, 0.75) - quantile(resid, 0.25))/1.34);
            f = resid./h;
            f = exp(-0.5 * f.^2) ;
            f = f./ sqrt(2*pi);
            f = f./h;
            
            fxxinv = bsxfun(@times,sqrt(f),xx);
            fxxt = fxxinv';
            fxxinv = inv(fxxt * fxxinv);
            
            cov = tau .* (1 - tau) .* fxxinv * xx' * xx * fxxinv;
        case 'wald'
            xxt = xx';
            xxi = xxt * xx;
            pz = sum(abs(resid) < sqrt(eps));
            h = max(p + 1, ceil(m * band(tau,m,inp.Results.hs)));
            ir = (pz + 1):(h + pz + 1);
            [~,Id]=sort(abs(resid));
            Id = resid(Id);
            ord_resid = sort(Id(ir));
            
            xt = ir'./(m - p);
            % sparsity = myqr3(xt', ord_resid, tau);
            sparsity = rq_fnm([ones(h+1,1) xt], ord_resid, 0.5);
            sparsity = sparsity(2);
            cov = sparsity.^2 .* tau .* (1 - tau).* inv(xxi) ;
    end
    serr = sqrt(diag(cov));
    tvalue = -abs(estimate./serr);
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
    pvalue = 2 * pt';
end
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


