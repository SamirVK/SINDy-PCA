clear all; close all; clc

%% Load in movies
load mri
s11 = load('cam1_1.mat');
s12 = load('cam1_2.mat');
s21 = load('cam2_1.mat');
s22 = load('cam2_2.mat');
s31 = load('cam3_1.mat');
s32 = load('cam3_2.mat');

%% Store the matrix struct fields 
% these have shape MxNx3xF
d11 = s11.vidFrames1_1;
d12 = s12.vidFrames1_2;
d21 = s21.vidFrames2_1;
d22 = s22.vidFrames2_2;
d31 = s31.vidFrames3_1;
d32 = s32.vidFrames3_2;

%% Play movies
mov = immovie(d32);
implay(mov)