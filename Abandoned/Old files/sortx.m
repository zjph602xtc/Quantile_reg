function [sortedx] = sortx( x,startx )
% input x is a table; start x is a string
% sortedx does NOT include startx

%%%%
indexall=x.Properties.VariableNames;
xn=length(indexall);
sortedx=cell(1,xn-1);

for i=1:xn-1
    xmain=x{:,startx};
    indexall=setdiff(indexall,startx);
    x=x(:,indexall);
    dis=cellfun(@disfun,num2cell(x{:,:},1));
    [~,I]=min(dis);
    sortedx(i)=indexall(I);
    startx=indexall{I};
end
end


function [dis] = disfun (x)
    xmain=evalin('caller','xmain');
    dis=sum(abs(x-xmain));
%     dis=sum((x-xmain).^2);
end

