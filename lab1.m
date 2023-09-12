num = [1280 640];
den = [1 24.2 1604.81 320.24 16];
sys = tf(num,den);

w= logspace(-4,3,100);
bode(sys,w)

