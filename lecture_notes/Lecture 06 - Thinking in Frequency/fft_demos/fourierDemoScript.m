% fourierDemoScript
close all
im_monaco = im2double(imread('../figs/im_monaco.png'));
im_town = im2double(imread('../figs/im_town.png'));
im_outdoor = im2double(imread('../figs/im_outdoor.png'));
im_beach = im2double(rgb2gray(imread('../figs/im_beach.jpg')));
im_straw = im2double(imread('../figs/straw.png'));

fourierDemo(im_monaco, false); pause;  clf;
fourierDemo(im_outdoor, false); pause;   clf;
fourierDemo(im_town, true); pause;  clf;
fourierDemo(im_beach, true); pause;   clf;


filterFourierDemo(im_town);

%filterFourierDemo(im_straw, 'straw');