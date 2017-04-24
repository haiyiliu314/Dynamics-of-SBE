f = importdata('f.dat');
psquared = importdata('psquared.dat');
dt = importdata('dt.dat');
tgrid = importdata('tgrid.dat');
time = dt * (1:(tgrid+1)) * 1e12;
plot(time,f)
hold on
plot(time,psquared)
xlabel('Time(ps)')
ylabel('Excitation')