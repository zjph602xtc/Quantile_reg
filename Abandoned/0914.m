n = 500;
p = 10;
x = normrnd(0,1,n,1);

z = normrnd(0,1,n,p);

bet0 = @(tau, alp00 , alp01, x, k){
        alp00 ./ p .* k + (alp01 + 0.5 .* x) .* log(tau);
    }

y = [];
for i = 1:n
    u = unifrnd(0,1,1,1);
    betas = [];
    for k = 1:p
        betas = [betas, cell2mat(bet0(u, -3.3,  0.86, x(i), k))];
    end
    y = [y; z(i,:)*betas'];
end

x = sort(x);

xx=[];
for j = (p + 2):(n - p - 2)
    xx = [xx 1 * (x <= x(j))];
end