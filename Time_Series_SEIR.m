%SEIR model, using Matlab to check Python
function Time_Series_SEIR()

%syms S E I R;

y0 = [90, 5, 5, 0];
tspan = [0 20];

[t, y] = ode45(@(t,y) model(t, y), tspan, y0);


figure;
    plot(t, y, 'LineWidth', 2);
    xlabel('Time (Days)');
    ylabel('Number of Individuals');
    legend('Susceptible (S)', 'Exposed (E)', 'Infectious (I)', 'Recovered (R)');
    title('SEIR Timeseries MATLAB');
    grid on;
end
 
function[dydt] = model(t, y)
S = y(1);
E = y(2);
I = y(3);
R = y(4);

gamma = 1;
%alpha = 0.29
alpha = 2.347510026817774; %transition time from low to high density; not infectivity!
beta = 3.374; %calculated from root-finding the imaginary values in the Jacobian
mu = 1/3e4;
Lambda = 1;

dSdt = Lambda - beta*S*I - mu*S;
dEdt = beta*S*I - alpha*E - mu*E;
dIdt = alpha*E - gamma*I - mu*I;
dRdt = gamma*I - mu*R;
dydt = [dSdt; dEdt; dIdt; dRdt];

end


