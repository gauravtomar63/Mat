Start=1;
End=1001;
%tr=100;
wl=50;
for tr=1:1001
for Start=1:End-wl
    
data3=Data(Start:Start+wl,tr);
ref3=Ref(Start:Start+wl,tr);

L=length(ref3);
win=hanning(L);
data4=data3(:).*win(:);
ref3=ref3(:).*win(:);


data1_fft=fft(data3);
ref_fft=fft(ref3);

cps1=ref_fft.*conj(data1_fft);
cps=cps1./abs(cps1);





%figure; plot(diff); ylim([-2.5 2]);
cc=xcorr(ref3,data3);
max_corr=max(cc);


for i=1:size(cc,1)
    if cc(i,1)==max_corr
        J=i;
    end
    
end

cc1=max_corr;
cc0=cc(J-1);
cc2=cc(J+1);

y0=cc0; y1=cc1; y2=cc2;
x0=J-1; x1=J; x2=J+1;
Y =[y0; y1; y2;];

X=[x0^2 x0 1; x1^2 x1 1; x2^2 x2 1];
Xinv=inv(X);

A=Xinv*Y;

a=A(1,:);
b=A(2,:);
c=A(3,:);
%y=a.*x^2+b.*x+c;

d=-b./(2.*a);
n=d-L;
timeshift(Start,tr)=n.*4;


qr=real(cps);
qi=imag(cps);
% result=0;
% for i=1:1001
% alpha=atan(qi(i,:)./qr(i,:));
% %thita=x.*n;
% syms = x;
% result(i,:) = int(x.*sin(alpha-x.*n),-pi,pi);
% end

%% gradient descent
ws=2*pi/L;
wnorm=-pi:ws:pi;
wnorm=wnorm(1:L);
qr1=sum(qr(:));
qi1=sum(qi(:));
%xk1=input('x10= ');
xk=n;

xk1=xk-.1;

%syms x1 alpha

% fw=0;
% for wx=0:2*pi
%     
%     fw=fw+2*wx.*(qr1.*sin(wx.*x1)-qi1.*cos(wx.*x1));
%     
% end    
p=0;
k=0;
%mk0=n;
xk1=xk-0.1;
d=.01;
deltam=0;
deltag=0;
gk1=1;
k_max=10;
while k<k_max
    fw=0;
for wx=1:L
    
    fw=fw+2*wx.*(qr(wx).*sin(wnorm(wx).*xk)-qi(wx).*cos(wnorm(wx).*xk));
    
end 
    %fwxk=subs(fw,x1,xk(1));
    
    if fw==0
        break;
    end
    
    deltam=xk-xk1;
    deltag=fw-gk1;
    d=deltam./deltag;
    xk_1=xk-d.*fw;
    %fprintf('x%d= %f \n',k, xk_1(1));
    if abs(xk-xk1)<.0001
        break;
    else
        xk1=xk;
        xk=xk_1;
        gk1=fw;
        k=k+1;
    end
    
end 
relativeshhft(Start,tr)=xk;
end
end
