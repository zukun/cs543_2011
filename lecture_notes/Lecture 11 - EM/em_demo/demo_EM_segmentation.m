function demo_EM_segmentation(im, K)

x = im(:);

[mu, sigma, prior, pm] = EM_gaussian(x, K);
xrange = 0:0.01:1;
figure(1), hold off
colors = 'rgbyck';
leg = cell(1, K);
for k = 1:K
    plot(xrange, prior(k)*normpdf(xrange, mu(k), sigma(k))*0.5, colors(mod(k-1,6)+1), 'linewidth', 2), 
    hold on;
    leg{k} = sprintf('{pi=%.2f, mu=%.2f, sigma=%.2f}', prior(k), mu(k), sigma(k));
end
legend(leg)

figure(2), hold off,
for k = 1:K
    pim = reshape(pm(:, k), size(im));
    subplot(1, K, k), imagesc(pim, [0 1]), axis image, axis off, colormap gray;
end

figure(3), hold off,
[spm, mi] = sort(pm, 2, 'descend');
for k = 1:K
    pim = zeros(size(im));
    pim((mi(:,1)==k) & (spm(:,1)>spm(:,2)+0.05)) = 1;
    subplot(1, K, k), imagesc(pim, [0 1]), axis image, axis off, colormap gray;
end

function [mu, sigma, prior, pm] = EM_gaussian(x, K)

x = x(:);
N = numel(x);

mu = zeros(K, 1);
sigma = zeros(K, 1);
minx = min(x); maxx = max(x);
for k = 1:K
    mu(k) = (0.1+0.8*rand(1))*(maxx-minx) + minx;
    sigma(k) = (rand(1)*0.9+0.1)*std(x);
end

% sx = sort(x);
% step = ceil(numel(x)/K);
% mu = zeros(K, 1);
% mu(:) = sx(ceil(step/2):step:end);
% sigma = zeros(K, 1);
% sigma(:) = rand(1)*std(x);
prior = zeros(K, 1);
prior(:) = 1/K;

pm = 1/K*ones(N, K);
oldpm = zeros(N, K);
logp = [];
maxiter = 200;
niter = 0;
while (mean(abs(pm(:)-oldpm(:)))>0.001) && (niter < maxiter)
   niter = niter+1;
    
    oldpm = pm;
    
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

    pdata = 0;
    for k = 1:K
        pdata = pdata + prior(k)*normpdf(x, mu(k), sigma(k));
    end
    logp(end+1) = mean(log(pdata));
    figure(4), hold off, plot(logp)
    legend(sprintf('Mean Log P(data) = %.3f', logp(end)))
    
end

