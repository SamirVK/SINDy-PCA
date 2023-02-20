clear all; close all; clc

%% Load data

% Load movies from file
s11 = load('cam1_1.mat');
s12 = load('cam1_2.mat');
s21 = load('cam2_1.mat');
s22 = load('cam2_2.mat');
s31 = load('cam3_1.mat');
s32 = load('cam3_2.mat');

% Load xpos and ypos arrays from file
s11_xpos = load('d11_xpos.mat');
s11_ypos = load('d11_ypos.mat');
s12_xpos = load('d12_xpos.mat');
s12_ypos = load('d12_ypos.mat');
s21_xpos = load('d21_xpos.mat');
s21_ypos = load('d21_ypos.mat');
s22_xpos = load('d22_xpos.mat');
s22_ypos = load('d22_ypos.mat');
s31_xpos = load('d31_xpos.mat');
s31_ypos = load('d31_ypos.mat');
s32_xpos = load('d32_xpos.mat');
s32_ypos = load('d32_ypos.mat');

% Store the matrix struct fields 
% These have shape MxNx3xF
d11 = s11.vidFrames1_1;
d12 = s12.vidFrames1_2;
d21 = s21.vidFrames2_1;
d22 = s22.vidFrames2_2;
d31 = s31.vidFrames3_1;
d32 = s32.vidFrames3_2;
% Only need up to minimum number of frames
frame_sizes = [size(d11,4); size(d12,4); size(d21,4); size(d22,4); 
    size(d31,4); size(d32,4)];
frames = min(frame_sizes);

% Store the array struct fields
d11_xpos = s11_xpos.xpos;
d11_ypos = s11_ypos.ypos;
d12_xpos = s12_xpos.xpos;
d12_ypos = s12_ypos.ypos;
d21_xpos = s21_xpos.xpos;
d21_ypos = s21_ypos.ypos;
d22_xpos = s22_xpos.xpos;
d22_ypos = s22_ypos.ypos;
d31_xpos = s31_xpos.xpos;
d31_ypos = s31_ypos.ypos;
d32_xpos = s32_xpos.xpos;
d32_ypos = s32_ypos.ypos;

%% Algorithm 1: Locate the yellow paint! 
% Using a yellow frame, we parse through each of the 6 videos looking 
% for the matrix that minimizes the Frobenius norm of the difference.

% Make a sort-of-pale yellow matrix
yellowArray = [255 255 255 255; 255 255 255 255; 255 255 255 255];
yellowArray(:,:,2) = [255 255 255 255; 255 255 255 255; 255 255 255 255];
yellowArray(:,:,3) = [100 100 100 100; 100 100 100 100; 100 100 100 100];
xpos = zeros(1,frames);
ypos = zeros(1,frames);

for k = 1:frames
    % Initialize the lowest score to a big number
    current_lowest_score = 1000000;
    winning_x_index = 0;
    winning_y_index = 0;
    for i = 1:478
        for j = 1:637
            % Switch d-- to desired video <<<<<<<<<<
            candidate_matrix = double(d32(i:i+2,j:j+3,:,k));
            candidate_score = norm(candidate_matrix-yellowArray,'Fro');
            if candidate_score <= current_lowest_score
                % Update tournament winners
                current_lowest_score = candidate_score;
                winning_x_index = j + 3;
                winning_y_index = i + 2;
            end
        end
    end
    % The gold medal goes to:
    xpos(k) = winning_x_index;
    ypos(k) = winning_y_index;
end
% Save the position vectors in 'd--_-pos.mat' <<<<<<<<<<
save('d32_xpos.mat','xpos');
save('d32_ypos.mat','ypos');

%% Play movies
load mri
mov = immovie(d22);
implay(mov)

%% Play video frames with centre position as blue dot
figure
for f = 1:frames
    imshow(d21(:,:,:,f))
    hold on
    plot(d21_xpos(f),d21_ypos(f),'.', 'MarkerSize',10)
    hold off
    pause(0.5);
end
