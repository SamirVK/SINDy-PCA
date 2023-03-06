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

%% Make the X-matrices (ideal and shaky) and perform PCA
X_ideal = [d11_xpos(1:frames); d11_ypos(1:frames); d21_xpos(1:frames); 
    d21_ypos(1:frames); d31_xpos(1:frames); d31_ypos(1:frames)];
X_shaky = [d12_xpos(1:frames); d12_ypos(1:frames); d22_xpos(1:frames); 
    d22_ypos(1:frames); d32_xpos(1:frames); d32_ypos(1:frames)];

% Subtract the mean from each row to obtain 0 mean
d11_xpos_av = mean(d11_xpos(1:frames));
d11_ypos_av = mean(d11_ypos(1:frames));
d21_xpos_av = mean(d21_xpos(1:frames));
d21_ypos_av = mean(d21_ypos(1:frames));
d31_xpos_av = mean(d31_xpos(1:frames));
d31_ypos_av = mean(d31_ypos(1:frames));
av_ideal = [d11_xpos_av; d11_ypos_av; d21_xpos_av; 
    d21_ypos_av; d31_xpos_av; d31_ypos_av];
X_ideal = (X_ideal - av_ideal);

d12_xpos_av = mean(d12_xpos(1:frames));
d12_ypos_av = mean(d12_ypos(1:frames));
d22_xpos_av = mean(d22_xpos(1:frames));
d22_ypos_av = mean(d22_ypos(1:frames));
d32_xpos_av = mean(d32_xpos(1:frames));
d32_ypos_av = mean(d32_ypos(1:frames));
av_shaky = [d12_xpos_av; d12_ypos_av; d22_xpos_av; 
    d22_ypos_av; d32_xpos_av; d32_ypos_av];
X_shaky = (X_shaky - av_shaky);

% Compute a rank reduced SVD decomposition of X_ideal and X_shaky
[U,S,V] = svd(X_ideal,'econ');
[Uu,Ss,Vv] = svd(X_shaky,'econ');

X_ideal_r2 = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
X_shaky_r2 = Uu(:,1:2)*Ss(1:2,1:2)*Vv(:,1:2)';

%% Apply the SINDy method to find coefficients of a linear system of ODEs.

fps = 20; % [1/s] "frames per second"
T = 1/fps; % [s] "seconds per frame (sampling period)"
L = (frames-51) * T; % [s] "length of video (sample length)"

t = 0:T:L; % time data
t = t(1,2:end);

% PCA using low-rank SVD
X = U(:,1:2)'*X_ideal;
Y = (X(:,2:end-50)-X(:,1:end-51))/T;
X = X(:,1:end-51);

% Populate the theta matrix by adding dictionary elements row-wise
theta = ones(1,frames-51); % constant term
% second order monomial dictionary:
theta(2:3,:) = X; % linear terms u and v
theta(4,:) = X(1,:).^2; % u^2
theta(5,:) = X(2,:).^2; % v^2
theta(6,:) = X(1,:).*X(2,:); % uv
theta(7,:) = X(1,:).^3; % u^3
theta(8,:) = X(2,:).^3; % v^3
theta(9,:) = X(1,:).*X(2,:).^2; % uv^2
theta(10,:) = X(1,:).^2.*X(2,:); % u^2v

Xi = Y*pinv(theta);

% Actual SINDy algorithm (from Jason Bramburger's SINDy.m)
% Sparsity parameter
lam = 0.5;

k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Loop over both rows 
    for ind = 1:2 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        Xi_new(ind,biginds) = Y(ind,:)*pinv(theta(biginds,:)); 
    end
    
    k = k + 1;
end

%% Print model (from Jason Bramburger's SINDy.m)

% Clear command window
clc

% Monomials up to degree 2
mons2 = ["" "u" "v" "u^2" "v^2" "uv"];
XiString = string(abs(Xi_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print du/dt model:
fprintf('du/dt = ')
bigcoeffs = abs(Xi_new(1,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(1,jnd),mons2(jnd)));
      else
        fprintf('+ %s', strcat(XiString(1,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dv/dt model:
fprintf('dv/dt = ')
bigcoeffs = abs(Xi_new(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(2,jnd) < 0  
        fprintf('- %s ', strcat(XiString(2,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(2,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% Try with integrals (from Jason Bramburger's SINDY.m)

fps = 20; % [1/s] "frames per second"
T = 1/fps; % [s] "seconds per frame (sampling period)"
L = (frames-51) * T; % [s] "length of video (sample length)"

t = 0:T:L; % time data
t = t(1,2:end);

X = U(:,1:2)'*X_ideal_r2; % Project X onto low-rank basis
X = X(:,1:end-51); % Truncate for better results
Y = X - X(:,1).*ones(size(X)); % Y is the X_tilde matrix

ThetaInt = zeros(size(theta)); % Initialize ThetaInt with all zeros

% Populate ThetaInt
ThetaInt(:,1) = T*theta(:,1);

for  n = 2:length(theta(1,:))
    ThetaInt(:,n) = ThetaInt(:,n-1) + T*theta(:,n);
end

% Solve the argmin problem
XiInt = Y*pinv(ThetaInt);

% Sparsity parameter
lam = 0.8;

k = 1;
XiInt_new = XiInt;
while sum(sum(abs(XiInt - XiInt_new))) > 0  || k == 1 
    
    XiInt = XiInt_new;
    
    % find small coefficients
    smallinds = (abs(XiInt) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    XiInt_new(smallinds) = 0;  
    
    % Loop over all 2 variables of the reduced model
    for ind = 1:2 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        XiInt_new(ind,biginds) = Y(ind,:)*pinv(ThetaInt(biginds,:)); 
    end
    
    k = k + 1;
end 

% Clear command window
clc

% Monomials up to degree 2
XiIntString = string(abs(XiInt_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print du/dt model:
fprintf('du/dt = ')
bigcoeffs = abs(XiInt_new(1,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if XiInt_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiIntString(1,jnd),mons2(jnd)));
      else
        fprintf('+ %s', strcat(XiIntString(1,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dv/dt model:
fprintf('dv/dt = ')
bigcoeffs = abs(XiInt_new(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if XiInt_new(2,jnd) < 0  
        fprintf('- %s ', strcat(XiIntString(2,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiIntString(2,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% Simulate the obtained model and compare with the data
t = 0:T:L; % time data
t = t(1,2:end);

figure(1)
plot(t,X(:,1:end))
axis([0 9 -120 120])
title('Paint can trajectory over time')
xlabel('Time (s)')
ylabel('Centred pixel data')
legend('Can position','Derivative')

t0 = t(1); 
tfinal = t(end);
uv0 = [42.1; 93.6];

[t,uv] = ode45(@NewtonHookeODE,[t0 tfinal],uv0);
[pks,locs] = findpeaks(uv(:,1));
dist = locs(2) - locs(1); %or whichever two peaks you want to know about
figure(2)
plot(t,uv)
axis([0 9 -140 140])
title('Simulated paint can trajectory over time')
xlabel('Time (s)')
ylabel('Centred pixel data')
legend('Can position','Derivative')


%% Algorithm 1: Locate the yellow paint! 
% Using a yellow frame, we parse through each of the 6 videos looking 
% for the matrix that minimizes the Frobenius norm of the difference.

% Make a pale yellow matrix
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
mov = immovie(d11);
implay(mov)

%% Play video frames with centre position as blue dot
figure
for f = 1:frames
    imshow(d22(:,:,:,f))
    hold on
    plot(d22_xpos(f),d22_ypos(f),'.', 'MarkerSize',10)
    hold off
    pause(0.5);
end

%% Functions

function duvdt = NewtonHookeODE(t,uv)

duvdt = [-3.0865-0.42511*uv(1)+3.1522*uv(2);
         -9.947-2.9197*uv(1)];
end
