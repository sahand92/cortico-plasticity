A_plus=1;
%A_minus=1; %CDP
A_minus=-1;%STDP
fmax=1000;
f = linspace(0,fmax,nw);
w=2*pi*f;

H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

tau_n=linspace(-10,0,200);
tau_p=linspace(0,10,200);
tau=linspace(-10,10,2*length(tau_n));
Htau=zeros(1,length(tau));
Htau(1:length(tau)/2)=-exp(tau_n);
Htau(length(tau)/2+1:length(tau))=exp(-tau_p);

Hw=(H0+1i*w.*tp*H1)./(1+w*tp).^2;
subplot(1,2,1)
plot(w,real(Hw))

hold on

plot(w,imag(Hw),'red')


hold on

plot(-w,-real(Hw))

hold on

plot(-w,-imag(Hw),'red')
xlabel('ω')
ylabel('')

hold on

legend('Re(H)','Im(Hw)')
xlabel('w (arbitrary units)')
ylabel('H (arbitrary units)')
subplot(1,2,2)

plot(tau,Htau)
xlabel('tau (arbitrary units)')
ylabel('H (arbitrary units)')
