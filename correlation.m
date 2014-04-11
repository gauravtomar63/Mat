N=160;
t=[0:N-1]'/N;
Y0=sin(2*pi*1*t); % reference signal with 0deg phase
Y45=sin(2*pi*1*t+(pi/4)); % signal with phase of 45 deg. 2 samples ahead with Y0
Y90=sin(2*pi*1*t+(pi/2)); % signal with phase of 90 deg. 4 samples a head with Y0
F(:,1)=(fft(Y0))/(N/2);
F(:,2)=(fft(Y45))/(N/2);
F(:,3)=(fft(Y90))/(N/2);
CS(:,1)=F(:,1).*conj(F(:,2));
CS(:,2)=F(:,1).*conj(F(:,3));
PEAK(:,1)=ifft(CS(:,1)); % look for corelation peak. It is on sample 3
PEAK(:,2)=ifft(CS(:,2));% look for corelation peak. It is on sample 5