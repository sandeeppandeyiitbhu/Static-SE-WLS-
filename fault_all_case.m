
%% IM Fault1
% [pf_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT1\Pflow_is_fault1_newmodel.csv');
% [p_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT1\Pinjection_is_fault1_newmodel.csv');
% [v_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT1\Vrms_is_fault1_newmodel.csv');
% time_data=p_data(:,1);
%% IM Fault2
% [pf_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT2\Pflow_is_fault2_newmodel.csv');
% [p_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT2\Pinjection_is_fault2_newmodel.csv');
% [v_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT2\Vrms_is_fault2_newmodel.csv');
% time_data=p_data(:,1);
%% IM Fault3
% [pf_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT3\Pflow_is_fault3_newmodel.csv');
% [p_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT3\Pinjection_is_fault3_newmodel.csv');
% [v_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT3\Vrms_is_fault3_newmodel.csv');
% time_data=p_data(:,1);
%% Grid connected Fault1
[pf_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT1\Pflow_gc_fault1_newmodel.csv');
[p_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT1\Pinjection_gc_fault1_newmodel.csv');
[v_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT1\Vrms_gc_fault1_newmodel.csv');
time_data=p_data(:,1);
%% Grid connected Fault2
% [pf_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT2\Pflow_gc_fault2_newmodel.csv');
% [p_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT2\Pinjection_gc_fault2_newmodel.csv');
% [v_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT2\Vrms_gc_fault2_newmodel.csv');
% time_data=p_data(:,1);
%% Grid connected Fault3
% [pf_data,~,~]=xlsread('C:\Users\pande\Desktop\New model\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT3\Pflow_gc_fault3_newmodel.csv');
% [p_data,~,~]=xlsread('C:\Users\pande\Desktop\New model\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT3\Pinjection_gc_fault3_newmodel.csv');
% [v_data,~,~]=xlsread('C:\Users\pande\Desktop\New model\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT3\Vrms_gc_fault3_newmodel.csv');
% time_data=p_data(:,1);
%% From P Injection data
PI3=p_data(:,7)/10;
PI5=p_data(:,9)/10;
PI7=p_data(:,11)/10;
QI3=p_data(:,19)/10;
QI5=p_data(:,21)/10;
QI7=p_data(:,23)/10;
%% From P Flow data
PF12=pf_data(:,2)/10;
PF23=pf_data(:,3)/10;
PF34=pf_data(:,4)/10;
PF37=pf_data(:,5)/10;
PF26=pf_data(:,6)/10;
PF65=pf_data(:,7)/10;
PF76=pf_data(:,8)/10;
QF12=pf_data(:,9)/10;
QF23=pf_data(:,10)/10;
QF34=pf_data(:,11)/10;
QF37=pf_data(:,12)/10;
QF26=pf_data(:,13)/10;
QF65=pf_data(:,14)/10;
QF76=pf_data(:,15)/10;
%% From V data
V1=v_data(:,4);
V1=V1/(13.2/sqrt(3));
V2=v_data(:,5);
V2=V2/(13.2/sqrt(3));
counter=0;
%% Actual Program Big Brain
for time_instant=time_data'
%%
   counter=counter+1;
 
           % Zdata
%         |Msnt |Type | Value | From | To | Rii | 
zdata7   = [ %---- Voltage Magnitude ------------%
            1     1     0.9887401391   1       0   9e-4;
            %-----------------------------------%
            %---- Real Power Injection ---------%
            2     2     0.970624852    3       0   1e-4;
            3     2     0.969349332    5       0   1e-4; 
            4     2     2.016405428    7       0   1e-4;
           %------------------------------------%
           %---- Reative Power Injection -------%
            5     3    0.430263378     3       0   1e-4;
            6     3   -0.932467604     5       0   1e-4; 
            7     3    1.347104622     7       0   1e-4;
           %------------------------------------%
           %------ Real Power Flow ------------- %
            8     4   -1.05086488      1       2   64e-6;
            9     4   -1.724532795     2       3   64e-6;
           10     4   -0.112682898     3       4   64e-6;
           11     4   -0.896995878     3       7   64e-6;
           12     4   -0.250640268     2       6   64e-6;
           13     4   -0.451520405     6       5   64e-6;
           14     4    0.991128501     7       6   64e-6;
           %------------------------------------%
           %------ Reactive Power Flow ---------%   
           15     5   0.87373507       1       2   64e-6;
           16     5   0.386070103      2       3   64e-6;
           17     5   0.579396034      3       4   64e-6;
           18     5  -0.663860084      3       7   64e-6;
           19     5   0.112710999      2       6   64e-6;
           20     5   0.660716265      6       5   64e-6;
           21     5   0.396254416      7       6   64e-6];
zdata7(2:end,3)=zdata7(2:end,3)/10;
%% 
%         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |
linedata7=[  1      2       0.07149   0.17840       0.0       1
             2      3       0.00596   0.01487       0.0       1
             2      6       0.01291   0.03221       0.0       1
             3      4       0.02284   0.05699       0.0       1
             3      7       0.00496   0.01239       0.0       1
             6      5       0.01688   0.04209       0.0       1
             6      7       0.01390   0.03469       0.0       1];
%% Creation of y bus
fb = linedata7(:,1);     % From bus number...
tb = linedata7(:,2);     % To bus number...
r = linedata7(:,3);      % Resistance, R...
x = linedata7(:,4);      % Reactance, X...
b = linedata7(:,5);      % Ground Admittance, B/2...
a = linedata7(:,6);      % Tap setting value..
z = r + i*x;            % Z matrix...
y = 1./z;               % To get inverse of each element...
b = i*b;                % Make B imaginary...

nbus = max(max(fb),max(tb));    % no. of buses...
nbranch = length(fb);           % no. of branches...
ybus = zeros(nbus,nbus);        % Initialise YBus...
 
 % Formation of the Off Diagonal Elements...
 for k=1:nbranch
     ybus(fb(k),tb(k)) = ybus(fb(k),tb(k))-y(k)/a(k);
     ybus(tb(k),fb(k)) = ybus(fb(k),tb(k));
 end
 
 % Formation of Diagonal Elements....
 for m =1:nbus
     for n =1:nbranch
         if fb(n) == m
             ybus(m,m) = ybus(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             ybus(m,m) = ybus(m,m) + y(n) + b(n);
         end
     end
 end
 %ybus;                  % Bus Admittance Matrix
 %zbus = inv(ybus);      % Bus Impedance Matrix

 %% wls
 num = nbus; 
%ybus = ybusppg(num); % Get YBus..nope we already have it
zdata = zdata7; % Get Measurement data..
    %%
fb = linedata7(:,1);
tb = linedata7(:,2);
b = linedata7(:,5);
nbus = max(max(fb),max(tb));    % no. of buses...
nbranch = length(fb);           % no. of branches...
bbus = zeros(nbus,nbus);

 for k=1:nbranch
     bbus(fb(k),tb(k)) = b(k);
     bbus(tb(k),fb(k)) = bbus(fb(k),tb(k));
 end
%%
bpq = bbus; % Get B data..
%nbus = max(max(zdata(:,4)),max(zdata(:,5))); % Get number of buses..
type = zdata7(:,2); % Type of measurement, Vi - 1, Pi - 2, Qi - 3, Pij - 4, Qij - 5, Iij - 6..
z =  [V1(counter) PI3(counter) PI5(counter) PI7(counter) QI3(counter) QI5(counter) QI7(counter) PF12(counter) PF23(counter) PF34(counter) PF37(counter) PF26(counter) PF65(counter) PF76(counter) QF12(counter) QF23(counter) QF23(counter) QF37(counter) QF26(counter) QF65(counter) QF76(counter)]';   %zdata7(:,3); % Measuement values..
fbus = zdata7(:,4); % From bus..
tbus = zdata7(:,5); % To bus..
Ri = diag(zdata7(:,6)); % Measurement Error..
V = ones(nbus,1); % Initialize the bus voltages..
del = zeros(nbus,1); % Initialize the bus angles..
E = [del(2:end); V];   % State Vector..
G = real(ybus);
B = imag(ybus);

vi = find(type == 1); % Index of voltage magnitude measurements..
ppi = find(type == 2); % Index of real power injection measurements..
qi = find(type == 3); % Index of reactive power injection measurements..
pf = find(type == 4); % Index of real powerflow measurements..
qf = find(type == 5); % Index of reactive powerflow measurements..

nvi = length(vi); % Number of Voltage measurements..
npi = length(ppi); % Number of Real Power Injection measurements..
nqi = length(qi); % Number of Reactive Power Injection measurements..
npf = length(pf); % Number of Real Power Flow measurements..
nqf = length(qf); % Number of Reactive Power Flow measurements..

iter = 1;
tol = 5;

while(tol > 1e-4)
% while(tol > 1e-4 && iter<4)
% while(tol > 2e-4)
    
    %Measurement Function, h
    h1 = V(fbus(vi),1);
    h2 = zeros(npi,1);
    h3 = zeros(nqi,1);
    h4 = zeros(npf,1);
    h5 = zeros(nqf,1);
    
    for i = 1:npi
        m = fbus(ppi(i));
        for k = 1:nbus
            h2(i) = h2(i) + V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
        end
    end
    
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:nbus
            h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
        end
    end
    
    for i = 1:npf
        m = fbus(pf(i));
        n = tbus(pf(i));
        h4(i) = -V(m)^2*G(m,n) - V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
    end
    
    for i = 1:nqf
        m = fbus(qf(i));
        n = tbus(qf(i));
        h5(i) = -V(m)^2*(-B(m,n)+bpq(m,n)) - V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
    end
    
    h = [h1; h2; h3; h4; h5];
    
    % Residue..
    r = z - h;
    
    % Jacobian..
    % H11 - Derivative of V with respect to angles.. All Zeros
    H11 = zeros(nvi,nbus-1);

    % H12 - Derivative of V with respect to V.. 
    H12 = zeros(nvi,nbus);
    for k = 1:nvi
        for n = 1:nbus
            if n == k
                H12(k,n) = 1;
            end
        end
    end

    % H21 - Derivative of Real Power Injections with Angles..
    H21 = zeros(npi,nbus-1);
    for i = 1:npi
        m = fbus(ppi(i));
        for k = 1:(nbus-1)
            if k+1 == m
                for n = 1:nbus
                    H21(i,k) = H21(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                H21(i,k) = H21(i,k) - V(m)^2*B(m,m);
            else
                H21(i,k) = V(m)* V(k+1)*(G(m,k+1)*sin(del(m)-del(k+1)) - B(m,k+1)*cos(del(m)-del(k+1)));
            end
        end
    end
    
    % H22 - Derivative of Real Power Injections with V..
    H22 = zeros(npi,nbus);
    for i = 1:npi
        m = fbus(ppi(i));
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H22(i,k) = H22(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                H22(i,k) = H22(i,k) + V(m)*G(m,m);
            else
                H22(i,k) = V(m)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
            end
        end
    end
    
    % H31 - Derivative of Reactive Power Injections with Angles..
    H31 = zeros(nqi,nbus-1);
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:(nbus-1)
            if k+1 == m
                for n = 1:nbus
                    H31(i,k) = H31(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                H31(i,k) = H31(i,k) - V(m)^2*G(m,m);
            else
                H31(i,k) = V(m)* V(k+1)*(-G(m,k+1)*cos(del(m)-del(k+1)) - B(m,k+1)*sin(del(m)-del(k+1)));
            end
        end
    end
    
    % H32 - Derivative of Reactive Power Injections with V..
    H32 = zeros(nqi,nbus);
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H32(i,k) = H32(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                H32(i,k) = H32(i,k) - V(m)*B(m,m);
            else
                H32(i,k) = V(m)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end
    end
    
    % H41 - Derivative of Real Power Flows with Angles..
    H41 = zeros(npf,nbus-1);
    for i = 1:npf
        m = fbus(pf(i));
        n = tbus(pf(i));
        for k = 1:(nbus-1)
            if k+1 == m
                H41(i,k) = V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
            else 
                if k+1 == n
                H41(i,k) = -V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                else
                    H41(i,k) = 0;
                end
            end
        end
    end
    
    % H42 - Derivative of Real Power Flows with V..
    H42 = zeros(npf,nbus);
    for i = 1:npf
        m = fbus(pf(i));
        n = tbus(pf(i));
        for k = 1:nbus
            if k == m
                H42(i,k) = -V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n))) - 2*G(m,n)*V(m);
            else 
                if k == n
                H42(i,k) = -V(m)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
                else
                    H42(i,k) = 0;
                end
            end
        end
    end
    
    % H51 - Derivative of Reactive Power Flows with Angles..
    H51 = zeros(nqf,nbus-1);
    for i = 1:nqf
        m = fbus(qf(i));
        n = tbus(qf(i));
        for k = 1:(nbus-1)
            if k+1 == m
                H51(i,k) = -V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            else 
                if k+1 == n
                H51(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
                else
                    H51(i,k) = 0;
                end
            end
        end
    end
    
    % H52 - Derivative of Reactive Power Flows with V..
    H52 = zeros(nqf,nbus);
    for i = 1:nqf
        m = fbus(qf(i));
        n = tbus(qf(i));
        for k = 1:nbus
            if k == m
                H52(i,k) = -V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n))) - 2*V(m)*(-B(m,n)+ bpq(m,n));
            else 
                if k == n
                H52(i,k) = -V(m)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                else
                    H52(i,k) = 0;
                end
            end
        end
    end
    
    % Measurement Jacobian, H..
    H = [H11 H12; H21 H22; H31 H32; H41 H42; H51 H52];
    
   % Gain Matrix, Gm..
    Gm = H'*inv(Ri)*H;
    
    %Objective Function..
    J = sum(inv(Ri)*r.^2); 
    
    % State Vector..
    dE = inv(Gm)*(H'*inv(Ri)*r);
    E = E + dE;
    del(2:end) = E(1:nbus-1);
    V = E(nbus:end);
    iter = iter + 1;
    tol = max(abs(dE));
end

CvE = diag(inv(H'*inv(Ri)*H)); % Covariance matrix..

Del = 180/pi*del;
E2 = [V Del]; % Bus Voltages and angles..
%disp('-------- State Estimation ------------------');
%disp('--------------------------');
%disp('| Bus |    V   |  Angle  | ');
%disp('| No  |   pu   |  Degree | ');
%disp('--------------------------');
for m = 1:nbus
    fprintf('%4g', m); fprintf('  %8.4f', V(m)); fprintf('   %8.4f', Del(m)); fprintf('\n');
end
disp('---------------------------------------------');
%% Estimated data
final_result(counter,:)=[time_instant V(1) Del(1) V(2) Del(2) V(3) Del(3) V(4) Del(4) V(5) Del(5) V(6) Del(6) V(7) Del(7)];
estimated=final_result(:,[2 4 6 8 10 12 14]);
estimated=estimated*(0.21 + 0.74)
bus1_E=estimated(:,1);
bus2_E=estimated(:,2);
bus3_E=estimated(:,3);
bus4_E=estimated(:,4);
bus5_E=estimated(:,5);
bus6_E=estimated(:,6);
bus7_E=estimated(:,7);

%% IM Fault1
[measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT1\Vrms_is_fault1_newmodel.csv');
t_data=measured_V_data(:,1);
measured_data=(measured_V_data(:,4:9));
bus1=measured_data(:,1)/(13.2/sqrt(3));
bus2=measured_data(:,5)/(13.2/sqrt(3));
bus3=measured_data(:,6)/(13.2/sqrt(3));
bus4=measured_data(:,2)/(13.2/sqrt(3));
bus6=measured_data(:,3)/(13.2/sqrt(3));
bus7=measured_data(:,4)/(13.2/sqrt(3));

%% IM Fault2
% [measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT2\Vrms_is_fault2_newmodel.csv');
% t_data=measured_V_data(:,1);
% measured_data=(measured_V_data(:,4:9));
% bus1=measured_data(:,1)/(13.2/sqrt(3));
% bus2=measured_data(:,5)/(13.2/sqrt(3));
% bus3=measured_data(:,6)/(13.2/sqrt(3));
% bus4=measured_data(:,2)/(13.2/sqrt(3));
% bus6=measured_data(:,3)/(13.2/sqrt(3));
% bus7=measured_data(:,4)/(13.2/sqrt(3));

%% IM Fault3
% [measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\IS_FAULT\IS_FAULT3\Vrms_is_fault3_newmodel.csv');
% t_data=measured_V_data(:,1);
% measured_data=(measured_V_data(:,4:9));
% bus1=measured_data(:,1)/(13.2/sqrt(3));
% bus2=measured_data(:,5)/(13.2/sqrt(3));
% bus3=measured_data(:,6)/(13.2/sqrt(3));
% bus4=measured_data(:,2)/(13.2/sqrt(3));
% bus6=measured_data(:,3)/(13.2/sqrt(3));
% bus7=measured_data(:,4)/(13.2/sqrt(3));

%% Grid connected Fault1
% [measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT1\Vrms_gc_fault1_newmodel.csv');
% t_data=measured_V_data(:,1);
% measured_data=(measured_V_data(:,4:9));
% bus1=measured_data(:,1)/(13.2/sqrt(3));
% bus2=measured_data(:,5)/(13.2/sqrt(3));
% bus3=measured_data(:,6)/(13.2/sqrt(3));
% bus4=measured_data(:,2)/(13.2/sqrt(3));
% bus6=measured_data(:,3)/(13.2/sqrt(3));
% bus7=measured_data(:,4)/(13.2/sqrt(3));

%% Grid connected Fault2
% [measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\same data\Microgrid_1_dataset_generator_rms-20240517T112617Z-001\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT2\Vrms_gc_fault2_newmodel.csv');
% t_data=measured_V_data(:,1);
% measured_data=(measured_V_data(:,4:9));
% bus1=measured_data(:,1)/(13.2/sqrt(3));
% bus2=measured_data(:,5)/(13.2/sqrt(3));
% bus3=measured_data(:,6)/(13.2/sqrt(3));
% bus4=measured_data(:,2)/(13.2/sqrt(3));
% bus6=measured_data(:,3)/(13.2/sqrt(3));
% bus7=measured_data(:,4)/(13.2/sqrt(3));

%%    GC Fault3
% % [measured_V_data,~,~]=xlsread('C:\Users\pande\Desktop\New model\Microgrid_1_dataset_generator_rms\GC_FAULT\GC_FAULT3\Vrms_gc_fault3_newmodel.csv');
% % t_data=measured_V_data(:,1);
% % measured_data=(measured_V_data(:,4:9));
% % bus1=measured_data(:,1)/(13.2/sqrt(3));
% % bus2=measured_data(:,5)/(13.2/sqrt(3));
% % bus3=measured_data(:,6)/(13.2/sqrt(3));
% % bus4=measured_data(:,2)/(13.2/sqrt(3));
% % bus6=measured_data(:,3)/(13.2/sqrt(3));
% % bus7=measured_data(:,4)/(13.2/sqrt(3));

end
% time_data=time_data(1:559)
% plot(time_data,bus2_E,'r',time_data,bus3_E,'g',time_data,bus4_E,'b',time_data,bus6_E,'m',time_data,bus7_E,'c')

hold on
% plot(t_data,bus2,'r',t_data,bus3,'g',t_data,bus4,'b',t_data,bus6,'m',t_data,bus7,'c')
plot(t_data,bus2,t_data,bus3,t_data,bus4,t_data,bus6,t_data,bus7)

% plot(timedata, x1, 'r', timedata, x2, 'g', timedata, x3, 'b');