%%
data = 'w45.big.opt';
types = {'dia' 'gra' 'pot' 'tot'};
for fc = types
    f = fc{:};
    A = dlmread(sprintf('%s.%s.time',data,f),' ');
    A(A==0) = NaN;
    [p,i] = sort(A(:,1));
    A = A(i,:);        
    op.(f).all = A;
    op.(f).p = A(:,1);
    op.(f).max = max(A(:,2:end),[],2);
    op.(f).min = min(A(:,2:end),[],2);
    op.(f).mean = mean(A(:,2:end),2);
    [p,i] = sort(op.(f).p);
end


%%
a = op.gra;
hold off;
plot( a.p, a.min(1)./a.min , '*-'); hold on;
plot( a.p, a.p, '-','color',[0.6 0.6 0.6], 'linewidth', 2);