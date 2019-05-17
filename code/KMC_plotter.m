L_loc = 100;
np = 4;
L = L_loc*np;
hInit = zeros(1, L);
hFinal = zeros(1, L);

for i = 0:np-1
    htemp = load(strcat('hInit0',num2str(i),'.txt'));
    hInit((L_loc*i+1):(L_loc*i + L_loc)) = htemp;
    htemp = load(strcat('hFinal0',num2str(i),'.txt'));
    hFinal((L_loc*i+1):(L_loc*i + L_loc)) = htemp;
end

figure(3);
plot(hFinal/L);
hold on;
plot(hInit/L);
legend('Final', 'Init');
hold off;


figure(4);
plot(hInit - hFinal);
