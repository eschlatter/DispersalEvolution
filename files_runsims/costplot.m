function D=costplot(v,nmax,M,K,B)
%calculates direct and indirect costs, and total benefit, of a given
%displacement strategy (for nmax=0,1,2)

D = length(v);

figure
hold on
bar(v,'w')
plot(1:D,M,'.-',1:D,K,'.-',1:D,B,'.-')
legend('kernel','mortality','kin competition','total benefit')
xlabel('Distance')
xticks(1:D)
xticklabels(string(0:D-1))
title(sprintf('nmax=%g',nmax))
hold off

end