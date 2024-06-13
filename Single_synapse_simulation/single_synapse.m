%% firing of a single symnpse with u and x updating 
% check how to incorporate decay of post synaptic Vm
% synaptic parameters

SimLength = 500;   %simulation length
tau_decay = 1500;    %recovery of synaptic resources 
tau_facil = 200;   %recovery time of utilisation factor 
tau_m   =  15;
U = 0.2;            %Baseline utilisation factor 
x_init = 1; 
Input = zeros(1,SimLength);                      
tps = [50; 80; 110; 140; 170; 200; 230; 260; 430];
Input(tps) = 1; 

U_log = zeros(1, SimLength);
X_log = zeros(1, SimLength);
Spikelog = zeros(1, SimLength);
V_log = zeros(1, SimLength);

%starting parameters
u = U; x = x_init; v = 0;
for t = 1:SimLength

    %log current values of u and x
    U_log(t) = u; X_log(t) = x;
    efficacy = u.*x;
    % only calculate post when receiving input
    if Input(t) >= 1
        spike = Input(t)*efficacy;
        Spikelog(t)    = spike;
    else
        spike = 0;
    end 

    V_log(t) = v;
    v = v + (1./tau_m) .* (spike-v);

    %update values
    u = u + (((U - u)/tau_facil) + U.*(1-u).*Input(t));
    x = x + (((1-x)/tau_decay) - (u.*x.*Input(t)));
end

figure;
subplot(3,1,1)
plot(1:SimLength,U_log,'b')
ylabel('u/x','FontSize',9)
hold on
plot(1:SimLength,X_log,'r')
% legend('u', 'x' ,'Location','southeast')

subplot(3, 1, 2)
plot(1:SimLength,Input,'k')
ylabel('Pre Vm','FontSize',9)
ylim([-0.25 1.25])

subplot(3,1,3)
plot(1:SimLength,V_log,'k')
ylim([-0.005 0.035])
ylabel('Post Vm','FontSize',9)
xlabel('Time (ms)','FontSize',9)
figure = gcf;
figure.Position = [100, 100, 300, 300];