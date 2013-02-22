function demo_linefit(npts, outliers, noise, method)
% demo_linefit(npts, outliers, noise, method)
% 
% Syntax examples
%   npts = 10; outliers = 0; noise = 0.05; method = {'lsq'};
%   demo_linefit(npts, outliers, noise, method);
% 
%   npts = 100; outliers = 80; noise = 0.01; method = {'hough', 'ransac'};
%   demo_linefit(npts, outliers, noise, method);

gt_m = -1;
gt_b = 1;

minx = 0; maxx = 1; %(1-gt_b)/gt_m;

if ischar(method), method = {method}; end

[x, y] = generateLine(npts, gt_m, gt_b, outliers/npts, noise, minx, maxx);

leg = {};
nl = 1;
figure(1), hold off, plot(x, y, '*r'), hold on;
leg{nl} = 'ground truth';
plot([minx maxx], [minx maxx]*gt_m+gt_b, '--k', 'linewidth', 2);
if any(strcmpi(method, 'lsq'))
    nl = nl+1;
    [m, b] = lsqfit(x, y);
    figure(1), plot([minx maxx], [minx maxx]*m + b, 'g', 'linewidth', 3); 
    leg{nl} = sprintf('Least Squares: m=%.2f b=%0.2f', m, b);
end
if any(strcmpi(method, 'tlsq'))
    nl = nl+1;
    [m, b] = total_lsqfit(x, y);
    figure(1), plot([minx maxx], [minx maxx]*m + b, 'b', 'linewidth', 3);     
    leg{nl} = sprintf('Total Least Squares: m=%.2f b=%0.2f', m, b);
end
if any(strcmpi(method, 'rlsq'))
    nl = nl+1;
    [m, b] = robust_lsqfit(x, y);
    figure(1), plot([minx maxx], [minx maxx]*m + b, '--g', 'linewidth', 3);     
    leg{nl} = sprintf('Robust Least Squares: m=%.2f b=%0.2f', m, b);
end
if any(strcmpi(method, 'hough'))
    nl = nl+1;
    [m, b] = houghfit(x, y);
    figure(1), plot([minx maxx], [minx maxx]*m + b, '--b', 'linewidth', 3);     
    leg{nl} = sprintf('Hough Transform: m=%.2f b=%0.2f', m, b);   
end
if any(strcmpi(method, 'ransac'))
    nl = nl+1;
    [m, b] = ransacfit(x, y);
    figure(1), plot([minx maxx], [minx maxx]*m + b, '-c', 'linewidth', 3);     
    leg{nl} = sprintf('RANSAC: m=%.2f b=%0.2f', m, b);   
end
axis([minx maxx sort([minx maxx]*gt_m+gt_b)]);
axis equal;
legend('', leg{:}, 'location', 'northeast')

function [x y] = generateLine(npoints, m, b, fractOutliers, noise_coef, minx, maxx)

x = minx + (maxx-minx)*rand(npoints, 1);
y = m*x+b;

% add random normal error in perpendicular direction
err = noise_coef*(randn(npoints, 1));
x = x+err*cos(atan(-1/m));
y = y+err*sin(atan(-1/m));

% for outliers, replace y-values with random values
nout = round(fractOutliers*npoints);
y(1:nout) = rand(nout, 1)*(max(y)-min(y))+min(y);

function [m, b] = lsqfit(x, y)
% y = mx + b
% find line that best predicts y given x
% minimize sum_i (m*x_i + b - y_i).^2
A = [x(:) ones(numel(x), 1)];
b = y(:);
p = A\b;
m = p(1);
b = p(2);

function [m, b, err] = total_lsqfit(x, y)
% ax + by + c = 0
% distance to line for (a^2+b^2=1): dist_sq = (ax + by + c).^2
A = [x(:)-mean(x) y(:)-mean(y)];
[v, d] = eig(A'*A);
p = v(:, 1); % eigenvector corr. to smaller eigenvalue

% get a, b, c parameters
a = p(1);
b = p(2);
c = -(a*mean(x)+b*mean(y));
err = (a*x+b*y+c).^2;

% convert to slope-intercept (m, b)
m = -a/b;
b = -c/b; % note: this b is for slope-intercept now

function [m, b] = robust_lsqfit(x, y)
% y = mx + b
% find line that best predicts y given x
% minimize sum_i (m*x_i + b - y_i).^2
% iterative robust fit
p = [0 ; 0];
sigma = 0;
sigmas = [];
for k = 1:5
    p = fminunc(@(p)geterr(p, x, y, sigma), p);
    warning off;
    err = sqrt((y-p(1)*x-p(2)).^2);
    sigma = median(err)*1.5;
    sigmas = [sigmas sigma];
end
m = p(1);
b = p(2);
figure(2), plot(sigmas)

function err = geterr(p, x, y, sigma)
err = (y-p(1)*x-p(2)).^2;
err = sum(err ./ (err + sigma.^2));

function [m, b] = houghfit(x, y)
% y = mx + b
% x*cos(theta) + y*sin(theta) = r
% find line that best predicts y given x
% minimize sum_i (m*x_i + b - y_i).^2
thetas = (-pi+pi/50):(pi/100):pi;
costhetas = cos(thetas);  
sinthetas = sin(thetas);
minr = 0; stepr = 0.005; maxr = 1; 

% count hough votes
counts = zeros(numel(thetas), (maxr-minr)/stepr+1);
for k = 1:numel(x)
    r = x(k)*costhetas + y(k)*sinthetas;
    inrange = find(r >= minr & r <= maxr); % only count parameters within the range of r
    rnum = round((r(inrange)-minr)/stepr)+1;
    ind = sub2ind(size(counts), inrange, rnum);
    counts(ind) = counts(ind) + 1;
end

% smooth the bin counts
counts = imfilter(counts, fspecial('gaussian', 5, 0.75));

% get best theta, rho and show counts
[maxval, maxind] = max(counts(:));
[thetaind, rind] = ind2sub(size(counts), maxind);
theta = thetas(thetaind);
r = minr + stepr*(rind-1);
figure(2), hold off, imagesc(counts'), axis image, colormap gray
xlabel('theta'); ylabel('rho')
title(sprintf('theta = %.2f  ;  rho = %.2f', theta, r));

% convert to slope-intercept
b = r/sin(theta);
m = -cos(theta)/sin(theta);

function [m, b] = ransacfit(x, y)
% y = mx + b
N = 200;
thresh = 0.03;

m = 0;
b = 0;
bestcount = 0;

for k = 1:N
    rp = randperm(numel(x));
    tx = x(rp(1:2));
    ty = y(rp(1:2));
    m = (ty(2)-ty(1)) ./ (tx(2)-tx(1));
    b = ty(2)-m*tx(2);
    
    nin = sum(abs(y-m*x-b)<thresh);
    if nin > bestcount
        bestcount = nin;
        inliers = (abs(y - m*x - b) < thresh);
    end
end

[m, b] = lsqfit(x(inliers), y(inliers));
figure(2), hold off, plot(x, y, 'r*'), hold on, plot(x(inliers), y(inliers), 'g*')
    






