function demo_EM_mog

mu_gt = [0 1];
sigma_gt = [2 0.5];
x(1:1000, 1) = randn(1000, 1)*sigma_gt(1)+mu_gt(1);
x(1001:2000, 1) = randn(1000, 1)*sigma_gt(2)+mu_gt(2);

% display ground truth
xrange = -10:0.01:10;
figure(1), hold off, plot(xrange, normpdf(xrange, mu_gt(1), sigma_gt(1))*0.5, 'b', 'linewidth', 2), hold on;
plot(xrange, normpdf(xrange, mu_gt(2), sigma_gt(2))*0.5, 'g', 'linewidth', 2)
pdata = 0;
for k = 1:2
    pdata = pdata + 0.5*normpdf(x, mu_gt(k), sigma_gt(k));
end
disp(['Mean Log P(data) for Ground Truth: ', num2str(mean(log(pdata)))])


%pause; %% XXXXX
[mu, sigma, prior] = EM_gaussian(x, 2);

figure(1), hold off, plot(xrange, normpdf(xrange, mu_gt(1), sigma_gt(1))*0.5, 'b', 'linewidth', 2), hold on;
plot(xrange, normpdf(xrange, mu_gt(2), sigma_gt(2))*0.5, 'g', 'linewidth', 2)
%plot(xrange, normpdf(xrange, mu_gt(1), sigma_gt(1))*0.5+normpdf(xrange, mu_gt(2), sigma_gt(2))*0.5, 'r');

plot(xrange, normpdf(xrange, mu(1), sigma(1))*prior(1), '--c'), hold on;
plot(xrange, normpdf(xrange, mu(2), sigma(2))*prior(2), '--r')
%plot(xrange, normpdf(xrange, mu(1), sigma(1))*prior(1)+normpdf(xrange, mu(2), sigma(2))*prior(2), '--r');
% disp(num2str(prior'))
% disp(num2str(mu'))
% disp(num2str(sigma'))

function [mu, sigma, prior] = EM_gaussian(x, K)

x = x(:);
N = numel(x);

% Random Initialization
mu = zeros(K, 1);
sigma = zeros(K, 1);
minx = min(x); maxx = max(x);
for k = 1:K
    mu(k) = (0.1+0.8*rand(1))*(maxx-minx) + minx;
    sigma(k) = (rand(1)*0.9+0.1)*std(x);
end
prior = zeros(K, 1);
prior(:) = 1/K;

pm = 1/K*ones(N, K);
oldpm = zeros(N, K);
logp = [];
while any(abs(pm-oldpm)>0.001) % convergence test
    
    oldpm = pm;
  
    % display current estimates
    xrange = minx:0.01:maxx;    
    leg = {};
    figure(2), hold off,
    colors = 'bgrcyk';
    for k = 1:K
        plot(xrange, prior(k)*normpdf(xrange, mu(k), sigma(k))*0.5, colors(mod(k-1,6)+1), 'linewidth', 2), 
        hold on;
        leg{k} = sprintf('{pi=%.2f, mu=%.2f, sigma=%.2f}', prior(k), mu(k), sigma(k));
    end
    legend(leg);
%     figure(2), hold off, plot(xrange, normpdf(xrange, mu(1), sigma(1))*prior(1), '--b', 'linewidth', 2), hold on;
%     plot(xrange, normpdf(xrange, mu(2), sigma(2))*prior(2), '--g', 'linewidth', 2)
%     legend(sprintf('{pi=%.2f, mu=%.2f, sigma=%.2f}', prior(1), mu(1), sigma(1)), ...
%         sprintf('{pi=%.2f, mu=%.2f, sigma=%.2f}', prior(2), mu(2), sigma(2)));        
    
    % estimate probability that each data point belongs to each component
    for k = 1:K
        pm(:, k) = prior(k)*normpdf(x, mu(k), sigma(k));
    end
    pm = pm ./ repmat(sum(pm, 2), [1 K]);
    
    % compute maximum likelihood parameters for expected densities
    for k = 1:K
        prior(k) = sum(pm(:, k))/N;
        mu(k) = sum(pm(:, k).*x) / sum(pm(:, k));
        sigma(k) = sqrt( sum(pm(:, k).*(x - mu(k)).^2)/sum(pm(:, k)));
    end
    

    % display likelihoods
    pdata = 0;
    for k = 1:K
        pdata = pdata + prior(k)*normpdf(x, mu(k), sigma(k));
    end
    logp(end+1) = mean(log(pdata));
    figure(3), hold off, plot(logp)
    legend(sprintf('Mean Log P(data) = %.3f', logp(end)))
    
    %pause; %% XXXXX
end
