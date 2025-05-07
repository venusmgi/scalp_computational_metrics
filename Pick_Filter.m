function [filter] = Pick_Filter(filterType,sampleFreq)
% Function for a switch case for picking a filter
% Rachel J. Smith (Lopouratory 2019)

Nyq = sampleFreq./2;

switch filterType
    case 'beta'
lowestFreq = 13; %Hz
secPerCycle = 1/lowestFreq;
%filter order should be at least 2 times(secPerCycle*sampleFreq)
n = round(13*(secPerCycle*sampleFreq)); 
filterBounds =  [13; 30];
% transision width should be 10-25% of the upper and lower frequency bounds
transWidth = [0.10; 0.2]; %having different transition width for different frequencyies; not related to the filterBounds order
f = [0 (1-transWidth(2))*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth(1))*filterBounds(2) Nyq]./Nyq;
a = [0 0 1 1 0 0];
filterKernel = firls(n,f,a);


% Alpha filter
    case 'alpha'
lowestFreq = 8; %Hz
secPerCycle = 1/lowestFreq;
%2 cycles at the lowest frequency analyzed
n = round(16*(secPerCycle*sampleFreq)); 
filterBounds =  [8; 13];
% transision width should be 10-25% of the upper and lower frequency bounds
transWidth = 0.1; 
f = [0 (1-transWidth)*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth)*filterBounds(2) Nyq]./Nyq;
a = [0 0 1 1 0 0];
filterKernel = firls(n,f,a);



% Theta filter
    case 'theta'
lowestFreq = 4; %Hz
secPerCycle = 1/lowestFreq;
%2 cycles at the lowest frequency analyzed
n = round(8*(secPerCycle*sampleFreq)); %cannot go higher than 8 because it will introduce a spike in filter response
filterBounds = [4;8];
transWidth = [0.15; 0.25]; %having different transition width for different frequencyies; not related to the filterBounds order
f = [0 (1-transWidth(2))*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth(1))*filterBounds(2) Nyq]./Nyq;
a = [0 0 1 1 0 0];
filterKernel = firls(n,f,a);


% Broadband filter
    case 'broadband'
lowestFreq = 2;
secPerCycle = 1/lowestFreq;
n = round(3.5*(secPerCycle*sampleFreq));
f = [0 0.25 1 55 58 Nyq]./Nyq;
a = [0 0 1 1 0 0];
filterKernel = firls(n,f,a);



    case 'delta'
lowestFreq = 1;
secPerCycle = 1/lowestFreq;
%2 cycles at the lowest frequency analyzed
%n must be at least 2 to 5 times ( 1/lowestFreq ) * fs ;
n = round(2*(secPerCycle*sampleFreq));
filterBounds = [1;4];
transWidth = 0.1; %having different transition width, to see which one will work better
f = [0 (1-transWidth)*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth)*filterBounds(2) Nyq]./Nyq;
a = [0 0 1 1 0 0];
filterKernel = fir2(n,f,a);


    otherwise 
        error('You did not choose one of the correct frequency band filters, you silly girl!')
end
end

