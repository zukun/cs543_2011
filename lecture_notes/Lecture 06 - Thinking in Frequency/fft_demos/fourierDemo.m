function fourierDemo(im, do_pause)

fftim = fft2(im);
fftim = fftshift(fftim);
fftimpow = log(abs(fftim+eps));
fftimpow = imfilter(fftimpow, fspecial('gaussian', [3 3], 0.5));
figure(1), hold off, imagesc(im), axis off, colormap gray, axis image
title('intensity image')
if do_pause
    pause;
end

sv = sort(fftimpow(:));  
minv = sv(1); maxv = sv(end);
minv = sv(round(0.005*numel(sv)));  %maxv = sv(round(0.999*numel(sv)));
figure(2), hold off, imagesc(fftimpow, [minv maxv]), axis off, colormap jet, axis image
colorbar
title('log fft magnitude')

% figure(3), hold off, imagesc(abs(fftim)), axis off, colormap jet, axis image
% title('fft magnitude')
% colorbar