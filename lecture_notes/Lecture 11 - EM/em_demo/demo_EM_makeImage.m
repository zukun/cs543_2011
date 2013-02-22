function [im, labels] = demo_EM_makeImage(K)

if ~exist('K', 'var')
    K = 1;
end

if K==1
    fg_p   = 1; 
    fg_mu  = 0.6;
    fg_std = 0.1;

    bg_p   = 1;
    bg_mu  = 0.4;
    bg_std = 0.1;
end

if K==2
    fg_p   = [0.5 0.5]; 
    fg_mu  = [0.25 0.75];
    fg_std = [0.05 0.05];

    bg_p   = 1;
    bg_mu  = 0.5;
    bg_std = 0.1;
end

%labels = zeros(100,100);
[x, y] = meshgrid(1:100, 1:100);
labels = (x-50).^2+(y-50).^2 <= 40.^2;
%labels(25:75, 25:75) = 1;
fgpix = find(labels==1);
bgpix = find(labels==0);
nfg = numel(fgpix);
nbg = numel(bgpix);

im = zeros(100, 100);

xfg = sample_mog(nfg, fg_p, fg_mu, fg_std);
im(fgpix) = xfg;

xbg = sample_mog(nbg, bg_p, bg_mu, bg_std);
im(bgpix) = xbg;


function x = sample_mog(N, p, mu, sigma)
% Sample from a mixture of univariate gaussians

x = zeros(N, 1);

K = numel(p);
x_all = zeros(numel(x), K);
for k = 1:K
    x_all(:, k) = randn(N, 1)*sigma(k)+mu(k);
end

for n = 1:N
    % sample a component
    k = find(rand(1)<cumsum(p), 1);
    x(n) = x_all(n, k);
end
    
