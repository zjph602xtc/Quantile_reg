
load('sample.mat')
tau=0.2;
Z=arrayfun(@(i)sprintf('Z%d',i),1:20,'unif',0);
X=arrayfun(@(i)sprintf('X%d',i),1:290,'unif',0);
y=samplemat{:,'Y'};

tic
nn=0;
for i=1:100
    inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,X(i)});
    [~,n ] = myqr([samplemat{:,[Z X(i)]} inter],y,tau);
    nn=nn+n';
end
toc

%%% add var
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,Z},y,tau);
TB=full(TB);
tic
nn=0;estall=[];
for i=1:100
    inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,X(i)});
    [~,n ] = multiaddmyqr([samplemat{:,X(i)} inter],A,TB,IB,c,nvar,y);
    nn=nn+n;
%     i
end
toc

%%% alter var
inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,X(1)});
[ A,TB,IB,c,nvar ] = mypreqr([inter samplemat{:,[X(1) Z]}],y,tau);
TB=full(TB);

tic
nn=0; estall=[];
iiiii=0;
for i=[1:100]
    inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,X(i)});
    [~,n ] = multialtermyqr([inter samplemat{:,X(i)}],A,TB,IB,c,nvar,y);
    nn=nn+n;
    if n~=0
        iiiii=iiiii+1;
    end
end
toc

%%% alter var sorted
startx='X1';
sortedx = sortx( samplemat(:,X),startx );

inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,startx});
[ A,TB,IB,c,nvar ] = mypreqr([inter samplemat{:,[startx Z]}],y,tau);
TB=full(TB);

tic
nn=0; estall=[];
iiiii=0;
for i=[1:100]
    inter = bsxfun(@times, samplemat{:,Z}, samplemat{:,sortedx(i)});
    [~,n ] = multialtermyqr([inter samplemat{:,sortedx(i)}],A,TB,IB,c,nvar,y);
    nn=nn+n;
    if n~=0
        iiiii=iiiii+1;
    end
end
toc