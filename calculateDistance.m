function distance = calculateDistance(target, duration, threshold, evalute)
[CF, AF, PF, f, Fs, L] = fftAnalyze(target, duration, threshold);
distance = evalute(CF, AF, PF, f, Fs, L);
end