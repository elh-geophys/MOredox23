% EL
% Updated 10-02-2023

% Figure combination


clear;

im_a = imread('EffvsD_H04_25th.png');
im_b = imread('EffvsD_H04_50th.png');
im_c = imread('EffvsD_H04_75th.png');

im_d = imread('EffvsD_N21_25th.png');
im_e = imread('EffvsD_N21_50th.png');
im_f = imread('EffvsD_N21_75th.png');


figure;
t = tiledlayout(2,3);
nexttile
imshow(im_a)
nexttile
imshow(im_b)
nexttile
imshow(im_c)
nexttile
imshow(im_d)
nexttile
imshow(im_e)
nexttile
imshow(im_f)

t.TileSpacing = 'none';
t.Padding = 'tight';