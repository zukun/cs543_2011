function swaptFFT_phase_magnitude_demo(im1, im2)

%% Initialization
if size(im1, 3)==3
    im1 = rgb2gray(im1);
end
if size(im2, 3)==3
    im2 = rgb2gray(im2);
end
im1 = imresize(im1, size(im2));
if max(size(im1))>640
    im1 = imresize(im1, 640/max(size(im1)));
    im2 = imresize(im2, 640/max(size(im2)));
end

%% Compute FFT and decompose to magnitude and phase
im1_fft = fft2(im1);
im1_fft_mag = abs(im1_fft); 
im1_fft_phase = angle(im1_fft);

im2_fft = fft2(im2);
im2_fft_mag = abs(im2_fft); 
im2_fft_phase = angle(im2_fft);

%% Display original images and after swapping phase
figure(1), imshow(im1);
figure(2), imshow(im2);
figure(3), imshow(ifft2(im1_fft_mag.*cos(im2_fft_phase)+1i*im1_fft_mag.*sin(im2_fft_phase)));
figure(4), imshow(ifft2(im2_fft_mag.*cos(im1_fft_phase)+1i*im2_fft_mag.*sin(im1_fft_phase)));
