clf
load 'clarinet_finger1.mat'
plot(f(1:16385), 20*log10(abs(Z(1:16385))), 'b')
hold on
load 'altosax_finger1.mat'
plot(f(1:16385), 20*log10(abs(Z(1:16385))), 'g')
xlim([0 5000])
hold off
grid
xlabel('Frequency (Hz)');
ylabel('Impedance / Zc (dB)');
legend('Clarinet Fingering 1', 'Alto Sax Fingering 1');

