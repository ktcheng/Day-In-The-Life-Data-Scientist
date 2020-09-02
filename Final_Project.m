%% Problem I Part B

% Declares array with N values and array to hold mu values
countArray = [10 20 50 100 200 300 500 1000 2000 10000 20000 30000 60000];
storageArray = [];

% Filters out NaN values in data.txt
arrayA = importdata('data.txt');
farrayA = arrayA(~isnan(arrayA));

n = 1; % Indexing variable
total = 0; % Keeps track of sum of data points

for i = 1:length(countArray)    
    while n < (countArray(i) + 1)
        total = total + farrayA(n); % Adds to total
        n = n + 1; % Increments index
    end
    
    mu = total / countArray(i); % Calculates mu
    storageArray = [storageArray, mu]; % Appends mu to array
end

% Plots mu_N v. N values
scatter(countArray, storageArray, 'filled'); grid on;
title('Mean of First N Non-NaN Elements v. N'); xlabel('Values of N'); 
ylabel('Mean of First N Non-NaN Elements');
%% Problem I Part C

% Declares array with N values and array to hold A_N values
countArray = [10 20 50 100 200 300 500 1000 2000 10000 20000 30000 60000];
accArray = [];

% Filters out NaN values in data.txt
arrayA = importdata('data.txt');
farrayA = arrayA(~isnan(arrayA));

for i = 1:length(countArray)
    total = 0; % Resets total on each iteration
    for k = 1:95241 % K_avail = 95241
        % Recall storageArray from Problem I Part B
        total = total + (farrayA(k) - storageArray(i)) .^ 2;
    end
    
    acc = total / length(farrayA); % Calculates A_N
    accArray = [accArray, acc]; % Appends A_N to array
end

% Plots A_N v. N Values
scatter(countArray, accArray, 'filled'); grid on; 
title('Accuracy of First N Non-NaN Elements v. N'); xlabel('Values of N'); 
ylabel('Accuracy of First N Non-NaN Elements');
%% Problem II Part A, C

% Declares array with N values and array to hold M_n values
secondN = [1 2 3 4 5 10 20 30];
meanArray = ones(8, 10000);

% Sets parameters for superimposed Gaussian RV (Part C)
x1 = [-2:0.01:10];
mu = 4;

for i = 1:length(secondN)
    % Generates N by 10000 samples of a uniform RV from [2, 6]
    randX = 2 + 4 * rand(secondN(i), 10000);
    
    % Appends M_n values to array
    meanArray(i, :) = sum(randX, 1) ./ secondN(i);
    
    subplot(2, 4, i);
    
    %%% PDF Plots
    histogram(meanArray(i, :), 'BinWidth', (1 / (secondN(i) + 1)), ...
        'normalization', 'pdf');
    titleString = strcat({'PDF N = '}, num2str(secondN(i)));
    title(titleString); xlabel('Sampled Value'); grid on;
    ylabel('Probability Density');
    
    %%% PDF Superimposed Gaussian
    hold on;
    varianceM = 4 / (3 * secondN(i)); % Sets the variance
    
    p1 = (1 / (sqrt(2 * pi * varianceM)));
    p2 = exp((-1 / 2) .* (((x1 - mu) .^ 2) / (varianceM)));
    fPdf = p1 .* p2; % PDF of superimposed Gaussian RV
    
    plot(x1, fPdf, 'LineWidth', 2); hold off; legend('M_{n}', 'Gaussian');
    
    %%% CDF Plots
    histogram(meanArray(i, :), 'BinWidth', (1 / (secondN(i) + 1)), ...
        'normalization', 'cdf');
    titleString = strcat({'CDF N = '}, num2str(secondN(i)));
    title(titleString); xlabel('Sampled Value'); grid on;
    ylabel('Cumulative Distribution');
    
    %%% CDF Superimposed Gaussian
    hold on;
    fCdf = normcdf(x1, mu, varianceM); % Recall variance from above
    
    plot(x1, fCdf, 'LineWidth', 2); hold off; legend('M_{n}', 'Gaussian');
end

% PDF Subplot Title
sgtitle('PDFs of Mean');

% CDF Subplot Title
sgtitle('CDFs of Mean');
%% Problem II Part D

% Declares array with N values and array to hold M_n values
secondN = [1 2 3 4 5 10 20 30];
meanArray = ones(8, 10000);
randDieX = [1 2 2 3 4 4 5]; % Weighted array of die values

% Sets parameters for superimposed Gaussian RV (Part C)
x1 = [-2:0.01:10];
mu = 3;

for i = 1:length(secondN)
    % Creates array to hold sampled values
    randXArray = ones(secondN(i), 10000);
    
    % Generates N by 10000 samples of die values
    for k = 1:secondN(i)
        % Generates 1 by 10000 samples of die values
        randX = datasample(randDieX, 10000); 
        randXArray(k, :) = randX; % Appends sampled values to array
    end
    
    % Appends M_n values to array
    meanArray(i, :) = sum(randXArray, 1) ./ secondN(i);
    
    subplot(2, 4, i);
    
    %%% PDF Plots
    histogram(meanArray(i, :), 'BinWidth', (1 / (secondN(i) + 1)), ...
        'normalization', 'pdf');
    titleString = strcat({'PDF N = '}, num2str(secondN(i)));
    title(titleString); xlabel('Sampled Value'); grid on;
    ylabel('Probability Density'); 
    
    %%% PDF Superimposed Gaussian
    hold on;
    varianceDie = 12 / (7 * secondN(i)); % Sets the variance
    
    p1 = (1 / (sqrt(2 * pi * varianceDie)));
    p2 = exp((-1 / 2) .* (((x1 - mu) .^ 2) / (varianceDie)));
    fPdf = p1 .* p2; % PDF of superimposed Gaussian RV
    
    plot(x1, fPdf, 'LineWidth', 2); hold off; legend('M_{n}', 'Gaussian');
    
    %%% CDF Plots
    histogram(meanArray(i, :), 'BinWidth', (1 / (secondN(i) + 1)), ...
        'normalization', 'cdf');
    titleString = strcat({'CDF N = '}, num2str(secondN(i)));
    title(titleString); xlabel('Sampled Value'); grid on;
    ylabel('Cumulative Distribution');
    
    %%% CDF Superimposed Gaussian  
    hold on;
    fCdf = normcdf(x1, mu, varianceDie); % Recall variance from above
    
    plot(x1, fCdf, 'LineWidth', 2); hold off; legend('M_{n}', 'Gaussian');
end

% PDF Subplot Title
sgtitle('PDFs of Mean');

% CDF Subplot Title
sgtitle('CDFs of Mean');
%% Problem III Part B, D

% Sets parameters from problem
samples = importdata('data_2.txt'); % Extracts data
mean0 = transpose([9 10]); % Initializes mu_0
mean1 = transpose([6 7]); % Initializes mu_1
tmean0 = transpose(mean0); % Transpose matrix of mean0
tmean1 = transpose(mean1); % Transpose matrix of mean1

% Initializes sigma matrix from problem
sigmaMatrix(1, :) = [1.15 0.1];
sigmaMatrix(2, :) = [0.1 0.5];
sigInverse = inv(sigmaMatrix); % Inverse matrix of sigmaMatrix

% Initializes arrays to hold class data points
class0 = [];
class1 = [];

% Initializes inequality terms
a = (1 / 2) * (tmean0 * sigInverse * mean0 - tmean1 * sigInverse * mean1);
a2 = a + (log(0.95 / 0.05)); % The 'a' term for part (d)

b = (sigInverse) * (mean1 - mean0);
tb = transpose(b); % Transpose matrix of b

% Initializes counters to keep track of class size
counter0 = 1;
counter1 = 1;

for i = 1:length(samples)
    x = samples(i, :); % Takes in one coordinate set
    
    %%% Part B Linear Inequality
    y = dot(transpose(b), x) + a;
    
    %%% Part D Linear Inequality
    y = dot(transpose(b), x) + a2;
    
    if (y < 0) % Classify as class 0 (based off of my signs)
        class0(counter0, :) = x; % Adds coordinates to class 0
        counter0 = counter0 + 1; % Increments class 0 counter
    else % Classify as class 1
        class1(counter1, :) = x; % Adds coordinates to class 1
        counter1 = counter1 + 1; % Increments class 1 counter
    end
end

% Percentage of samples from class 0
disp(((counter0 - 1) / length(samples)) * 100);

scatX0 = class0(:, 1); % X-coordinates from class 0
scatY0 = class0(:, 2); % Y-coordinates from class 0
sct0 = scatter(scatX0, scatY0, 'red'); hold on; % Plots scatter of class 0

scatX1 = class1(:, 1); % X-coordinates from class 1
scatY1 = class1(:, 2); % Y-coordinates from class 1
sct1 = scatter(scatX1, scatY1, 'blue'); hold on; % Plots scatter of class 1

legend([sct0, sct1], {'Class 0', 'Class 1'}); % Plot legend

% Initializes xy values for contour equation
X = [-10:0.01:10];
Y = [-10:0.01:10];

%%% Part B Contour Equation
y1 = @(X, Y) -((tb(1) * X) + tb(2) * Y) - a;

%%% Part D Contour Equation
y1 = @(X, Y) -((tb(1) * X) + tb(2) * Y) - a2;

% Plots contour line
fcontour(y1, 'LevelList', 0); title('Linear Discriminant Analysis'); grid on;
xlabel('X Values (from Data\_2)'); ylabel('Y Values (from Data\_2)');
%% Problem III Part F

% Sets parameters from problem
samples = importdata('data_3.txt'); % Extracts data
mean0 = transpose([9 10]); % Initializes mu_0
mean1 = transpose([6 7]); % Initializes mu_1
tmean0 = transpose(mean0); % Transpose matrix of mean0
tmean1 = transpose(mean1); % Transpose matrix of mean1

% Initializes sigma matrices from problem
sigmaMatrix0(1, :) = [1.15 0.1]; 
sigmaMatrix0(2, :) = [0.1 0.5];

sigmaMatrix1(1, :) = [0.2 0.3];
sigmaMatrix1(2, :) = [0.3 2];

sigInverse0 = inv(sigmaMatrix0); % Inverse matrix of sigmaMatrix0
sigInverse1 = inv(sigmaMatrix1); % Inverse matrix of sigmaMatrix1

% Initializes arrays to hold class data points
class0 = [];
class1 = [];

% Initializes inequality terms
c = (1 / 2) * (inv(sigmaMatrix0) - inv(sigmaMatrix1));
a = (1 / 2) * (tmean0 * sigInverse0 * mean0 - tmean1 * sigInverse1 * mean1);
a2 = a + (1 / 2) * log((det(sigmaMatrix0)) / (det(sigmaMatrix1)));

b = (sigInverse1 * mean1 - sigInverse0 * mean0);
tb = transpose(b); % Transpose matrix of b

% Initializes counters to keep track of class size
counter0 = 1;
counter1 = 1;

for i = 1:length(samples)
    x = samples(i, :); % Takes in one coordinate set
    tx = transpose(x); % Transpose matrix of x
    y = dot(transpose(b), x) + a2; % Latter terms of quadratic inequality
    
    % First term of quadratic inequality
    cDot = (c * transpose(x));
    cDot2 = dot(transpose(x), cDot);
    y = y + cDot2;
    
    if (y < 0) % Classify as class 0 (based off of my signs)
        class0(counter0, :) = x; % Adds coordinates to class 0
        counter0 = counter0 + 1; % Increments class 0 counter
    else
        class1(counter1, :) = x; % Adds coordinates to class 1
        counter1 = counter1 + 1; % Increments class 1 counter
    end
end

% Percentage of samples from class 0
disp(((counter0 - 1) / length(samples)) * 100);

scatX0 = class0(:, 1); % X-coordinates from class 0
scatY0 = class0(:, 2); % Y-coordinates from class 0
sct0 = scatter(scatX0, scatY0, 'red'); hold on; % Plots scatter of class 0

scatX1 = class1(:, 1); % X-coordinates from class 1
scatY1 = class1(:, 2); % Y-coordinates from class 1
sct1 = scatter(scatX1, scatY1, 'blue'); hold on; % Plots scatter of class 1

legend([sct0, sct1], {'Class 0', 'Class 1'}); % Plot legend

% Initializes xy values for contour equation
X = [-10:0.01:10];
Y = [-10:0.01:10];

% Inputs xy values into a single matrix
xyMatrix(:, 1) = X;
xyMatrix(:, 2) = Y;

% Calculates quadratic term of quadratic inequality
a3 = c * transpose(xyMatrix);
a4 = dot(transpose(xyMatrix), a3);
c11 = c(1, 1); c12 = c(1, 2); c21 = c(2, 1); c22 = c(2, 2);

%%% Part F Contour Equation
y1 = @(X, Y) (-((tb(1) * X) + tb(2) * Y) - a2 - ...
    (c11 * (X .^ 2) + c12 * (Y .* X) + c21 * (X .* Y) + c22 * (Y .^ 2)));

% Plots contour line
fcontour(y1, 'LevelList', 0, 'LineWidth', 2); 
xlabel('X Values (from Data\_3)'); ylabel('Y Values (from Data\_3)');
title('Quadratic Discriminant Analysis'); grid on;