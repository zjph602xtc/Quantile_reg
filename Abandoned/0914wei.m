
tau=0.5;

tic
nn=0;
for i=1:7:437
    inter = bsxfun(@times, z, xx(:,i));
    [~,n ] = myqr([inter xx(:,i) z],y,tau);
    nn=nn+n';
    i
end
toc

%%% add var
[ A,TB,IB,c,nvar ] = mypreqr(z,y,tau);
TB=full(TB);
tic
nn=0;
for i=1:437
    inter = bsxfun(@times, z, xx(:,i));
    [~,n ] = multiaddmyqr([inter xx(:,i)],A,TB,IB,c,nvar,y);
    nn=nn+n;

end
toc


%%% alter var sorted
inter = bsxfun(@times, z, xx(:,1));
[ A,TB,IB,c,nvar ] = mypreqr([inter xx(:,1) z],y,tau);
TB=full(TB);

tic
nn=0; estall=[];
iiiii=0;
for i=2:400
    inter = bsxfun(@times, z, xx(:,i));
    [~,n, A,TB,IB] = multialtermyqr([inter xx(:,i)],A,TB,IB,c,nvar,y);
    nn=nn+n;
%     i
end
toc