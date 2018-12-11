%plotting errorbars

h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');

%coded vs uncoded
a1=y{1}; b1=y{2};

a2=y{6}; b2=y{5};

a3=y{4}; b3=y{3};

iset=[1,5,10,15];
errorbar(x{1}(iset),a3(iset),zeros(size(x{1}(iset))),a2(iset)-a3(iset),'.r','linewidth',1,'CapSize',12)

errorbar(x{1}(iset),a3(iset),zeros(size(x{1}(iset))),a1(iset)-a3(iset),'.r','linewidth',1,'CapSize',12)

iset=[1,5,10,15]+2;
errorbar(x{1}(iset),b3(iset),zeros(size(x{1}(iset))),b2(iset)-b3(iset),'.k','linewidth',1,'CapSize',12)

errorbar(x{1}(iset),b3(iset),zeros(size(x{1}(iset))),b1(iset)-b3(iset),'.k','linewidth',1,'CapSize',12)



figure
for i=1
plot(x{i},y{i})
end