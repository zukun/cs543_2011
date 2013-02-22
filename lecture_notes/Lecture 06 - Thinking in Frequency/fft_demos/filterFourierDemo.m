function filterFourierDemo(im, special)

if ~exist('special', 'var')
    special = [];
end

if strcmp(special, 'straw')
    fil = fspecial('gaussian', 31, 4);
    displayAll(im, fil, 'gaussian')
    fil = ones(9, 9)/81;  
    displayAll(im, fil, 'box')
else
    fil = fspecial('sobel')'; %[0 0 0 ; -1 0 1 ; 0 0 0];
    displayAll(im, fil, 'sobel')

    fil = fspecial('gaussian', 31, 4);
    displayAll(im, fil, 'gaussian')

    fil = zeros(size(fil)); fil(15, 15) = 1; fil = fil - fspecial('gaussian', 31, 4); 
    displayAll(im, fil, 'log')

    fil = ones(9, 9)/81;  
    displayAll(im, fil, 'box')
end
% fil = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];  
% displayAll(im, fil)


function displayAll(im, fil, fname)

fftim = fft2(im);
fftim = fftshift(fftim);
fftimpow = log(abs(fftim+eps));
figure(1), hold off, imagesc(im), axis off, colormap gray, axis image
title('intensity image')
sv = sort(fftimpow(:));  
minv = sv(1); maxv = sv(end);
minv = sv(round(0.005*numel(sv)));  %maxv = sv(round(0.999*numel(sv)));
figure(2), hold off, imagesc(fftimpow, [minv maxv]), axis off, colormap jet, axis image
title('log fft magnitude of image')

fil = padarray(fil, [1 1], 0, 'both');
figure(3), hold off, imagesc(fil), axis off, colormap gray, axis image
title(['filter: ' fname])
fil2 = zeros(size(im));
fil2((1:size(fil, 1))+ floor((size(im,1)-size(fil,1))/2), ...
    (1:size(fil, 2))+ floor((size(im,2)-size(fil,2))/2)) = fil;
fftfil = fft2(fil2);
fftfil = fftshift(fftfil);
fftfilpow = abs(fftfil+eps);
figure(4), hold off, imagesc(fftfilpow), axis off, colormap jet, axis image, colorbar
title(['filter: ' fname])

im2 = imfilter(im, fil, 'same');
fftim2 = fft2(im2);
fftim2 = fftshift(fftim2);
fftim2pow = log(abs(fftim2+eps));
figure(5), hold off, imagesc(im2), axis off, colormap gray, axis image
title(['filtered image'])
figure(6), hold off, imagesc(fftim2pow, [minv maxv]), axis off, colormap jet, axis image, colorbar
title(['log fft magnitude of filtered image'])
% fftim3 = fft2(im).*fft2(fil2);
% fftim3 = fftshift(fftim3);
% fftim3pow = log(abs(fftim3+eps));
% figure(7), hold off, imagesc(fftim3pow, [minv maxv]), axis off, colormap jet, axis image, colorbar
pause;


