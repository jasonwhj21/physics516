estimates = [3.146721,3.071978,3.156832,3.148029,3.142598,3.141363];
errors = [0.208166,0.06418303,0.02029343,0.006455769,0.002038160,0.0006436586];
trials = [10,100,1000,10000,100000,1000000];

subplot(1,2,1)
errorbar(trials,estimates,errors,"both","o")
set(gca,'xscale','log')
xlabel('Monte Carlo Trials') 
ylabel('Pi Estimate')

subplot(1,2,2)
stdv = [0.211797,0.065181,0.020283,0.006990,0.001965, 0.000611];
logstdv = log10(stdv);
logtrials = log10(trials);
xs = [0,6]
ys = -0.5064*xs -0.1657
scatter(logtrials, logstdv,50,'o','filled','k');
hold on
scatter(log10(trials),log10(errors),50,'o','filled','r')
plot(xs,ys)
legend('Actual stdv','Unbiased Estimate stdv','-0.5064*x -0.1657, p=-0.5064')
xlabel('log(numTrials)')
ylabel('Standard Deviation')