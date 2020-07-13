function [CF, AF, PF, f, Fs, L] = fftAnalyze(target, duration, threshold)
if nargin == 2
    threshold = 1e-3;
elseif nargin ~= 3
    error('2 or 3 arguments are required to perform fftAnalyze.');
end
L = size(target, 2);
Fs = 1/(duration/L);
Y = fft(target, L, 2);
P2 = Y/L;
P1 = P2(:, 1:(L/2 + 1));
P1(:, 2:end - 1) = 2*P1(:, 2:end - 1);
f = Fs*(0:(L/2))/L;
AF = abs(P1);
maxAF = max(AF, [], 2);
PF = angle(P1).*(AF > maxAF*threshold);
CF(:, :, 1) = -AF.*sin(PF);
CF(:, :, 2) = AF.*cos(PF);
end
