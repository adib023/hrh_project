clc
clear all

w = 1.7565;
N = 10;
final_t = 500;
Npts = 2^10;
x = zeros(Npts,1);
t = linspace(0,final_t,Npts);

for p = 1:Npts
    if t(p) > pi*N/w 
        break
    end
    x(p) =  sin(w*t(p))*sin(w/N*t(p));
end

y = fft(x);
p2 = abs(y/Npts);
p1 = p2(1:Npts/2+1);
dt = final_t/Npts;
f = 1/dt*(0:Npts/2)/Npts*2*3.1416;
figure
subplot(2,1,1)
plot(t,x);
xlabel("$t$","interpreter","latex","FontSize",16);
ylabel("$F(t)$","interpreter","latex","FontSize",16);
subplot(2,1,2)
plot(f,p1);
xlabel("$\omega$","interpreter","latex","FontSize",16);
ylabel("$F(\omega)$","interpreter","latex","FontSize",16);
