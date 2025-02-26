function [filter] = Pick_Filter(filterType,sampleFreq)
% Function for a switch case for picking a filter
% Rachel J. Smith (Lopouratory 2019)

Nyq = sampleFreq./2;

switch filterType
    case 'beta'
lowestFreq = 14; %Hz
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq)); %2 cycles at the lowest frequency analyzed
f = [0 11 14 30 35 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);


% Alpha filter
    case 'alpha'
lowestFreq = 8; %Hz
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq)); %2 cycles at the lowest frequency analyzed
f = [0 4 8 12 14 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);


% Theta filter
    case 'theta'
lowestFreq = 3; %Hz
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq)); %2 cycles at the lowest frequency analyzed
f = [0 3 4 7 9 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);


% Broadband filter
    case 'broadband'
lowestFreq = 2;
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq));
f = [0 0.25 1 55 58 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);


% No-Muscle Filter (filter for frequencies associated with muscle
% artifact)
    case 'muscle'
lowestFreq = 15;
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq));
f = [0 15 20 50 55 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);

    case 'delta'
lowestFreq = 1;
secPerCycle = 1/lowestFreq;
n = ceil(2*(secPerCycle*sampleFreq));
f = [0 0.2 1 4 6 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filter = firls(n,f,a);

% % If you want to check the filter:
% [h,w] = freqz(filter,1,500,2);
% 
% figure;
% plot(f,a,'r.-');
% hold on;
% plot(w,abs(h));

    otherwise 
        error('You did not choose one of the correct frequency band filters, you silly girl!')
end
end

