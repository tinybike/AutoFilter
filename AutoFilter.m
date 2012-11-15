% AutoFilter.m
% Estimates the autoregression kernel of a noisy input signal
% using a simple low-pass filter.

% Variable definitions:
% tests: number of different corrupted signals to analyze
% its: maximum filter power to test on each signal
% noisemax: number of different noise levels to test
% nstd: standard deviation of white noise
% frac_success: ratios of successful kernel estimates
% median_success_n: median filter powers for successful kernel estimates
% SNR: signal-to-noise ratio
tests = 100;
its = 20;
noisemax = 200;
nstd = linspace(0.01,3.5,noisemax);
frac_success = zeros(1,noisemax);
median_success_D = zeros(1,noisemax);
median_success_p = zeros(1,noisemax);
median_success_n = zeros(1,noisemax);
SNR = zeros(1,noisemax);

% Define signal (sinc function)
t = 0.01:0.5:50;
S = sin(t)./t;
T = length(t);

% Calculate AR kernel
dS = S(2:T) - S(1:T-1);
P = sparse(toeplitz(S(1:T-1),cat(2,S(1),zeros(1,T-2))));
K = P\dS';

% "Running average" low-pass filter
f = sparse(cat(2,[2 1],zeros(1,T-2)));
F = toeplitz(f)/4;

for k = 1:noisemax
    
    % Calculate SNR
    SNR(k) = 1/mean(nstd(k)*abs(randn(1,10000)));
    
    % Run tests to see how often the kernel is successfully calculated
    success = 0;
    clear success_n
    clear success_D
    clear success_p
    for j = 1:tests
        
        % Add white noise to signal
        E = nstd(k)*randn(1,T).*S;
        SE = S + E;
        dSE = SE(2:T) - SE(1:T-1);
        
        % Loop to optimize the filter power
        KS = zeros(its,3);
        KEf = zeros(T-1,its);
        accept = 0;
        clear n
        for i = 1:its
            
            % Filter signals and estimate kernel
            Fn = F^i;
            SEf = Fn*SE';
            dSEf = SEf(2:T) - SEf(1:T-1);
            dSEf = Fn(1:T-1,1:T-1)*dSEf;
            PEf = sparse(toeplitz(SEf(1:T-1),cat(2,SEf(1),zeros(1,T-2))));
            KEf(:,i) = PEf\dSEf;
            KEf(:,i) = Fn(1:T-1,1:T-1)*KEf(:,i);
            
            % Two-sample KS test of estimated vs. exact kernel
            [h,p,D] = kstest2(K,KEf(:,i));
            KS(i,:) = [h p D];
            
        end
        
        % Pick optimal filter power (using KS statistic)
        [n,n] = min(KS(:,3));
        success = success + ~KS(n,1);
        if ~KS(n,1)
            success_n(success) = n;
            success_p(success) = KS(n,2);
            success_D(success) = KS(n,3);
        end
        
    end
    
    frac_success(k) = success/tests;
    if success
        median_success_n(k) = median(success_n);
        median_success_p(k) = median(success_p);
        median_success_D(k) = median(success_D);
    end
    
end

SR = [1./SNR' frac_success'];
MFP = [1./SNR' median_success_n'];
MSP = [1./SNR' median_success_p'];
MSD = [1./SNR' median_success_D'];
SR = sortrows(SR,1);
MFP = sortrows(MFP,1);
MSP = sortrows(MSP,1);
MSD = sortrows(MSD,1);

figure
clf
plot(SR(:,1),SR(:,2))
xlabel('1/SNR')
ylabel('Success ratio (p > 0.05)')

figure
clf
plot(MFP(:,1),MFP(:,2))
xlabel('1/SNR')
ylabel('Median filter power')

figure
clf
plot(MSP(:,1),MSP(:,2))
xlabel('1/SNR')
ylabel('Median p-value')

figure
clf
plot(MSD(:,1),MSD(:,2))
xlabel('1/SNR')
ylabel('Median KS statistic')