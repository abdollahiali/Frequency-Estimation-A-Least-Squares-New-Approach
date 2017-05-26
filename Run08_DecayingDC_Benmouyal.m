clc
clear
fontSize = 10;
format long

t_sample = 0.13;
f0 = 50;
w0 = 2*pi*f0;

N = 20;
fs = N*f0;
dt = 1/fs;

Arms = 1;
Am = Arms*sqrt(2);

for iteration=1:2
    % begin parameters of decaying DC
    switch(iteration)
        case 1

            A_DC = 0.2*Arms;
            A_DC_a = 1*A_DC;
            A_DC_b = 0*A_DC;
            A_DC_c = 0*A_DC;
            T = 1/f0;
            coeff_tau = 1.5;
            coeff_tau_a = 1*coeff_tau;
            coeff_tau_b = 1*coeff_tau;
            coeff_tau_c = 1*coeff_tau;
            % coeff_tau_a = 0.5;
            % coeff_tau_b = 1;
            % coeff_tau_c = 3;

            tau_a = coeff_tau_a*T;
            r_a = exp(-dt/tau_a);

            tau_b = coeff_tau_b*T;
            r_b = exp(-dt/tau_b);

            tau_c = coeff_tau_c*T;
            r_c = exp(-dt/tau_c);
        case 2
            A_DC = 0.2*Arms;
            A_DC_a = 1*A_DC;
            A_DC_b = 1*A_DC;
            A_DC_c = 1*A_DC;
            T = 1/f0;
            coeff_tau = 1.5;
            coeff_tau_a = 1*coeff_tau;
            coeff_tau_b = 1*coeff_tau;
            coeff_tau_c = 1*coeff_tau;
            % coeff_tau_a = 0.5;
            % coeff_tau_b = 1;
            % coeff_tau_c = 3;

            tau_a = coeff_tau_a*T;
            r_a = exp(-dt/tau_a);

            tau_b = coeff_tau_b*T;
            r_b = exp(-dt/tau_b);

            tau_c = coeff_tau_c*T;
            r_c = exp(-dt/tau_c);
    end

    % end parameters of decaying DC


    tf = 5*coeff_tau*T-dt;

    nFilter = 0;
    f1 = 20;
    f2 = 90;

    h3Coeff = 0;% 0.05;
    h5Coeff = 0;%0.03;

    % do not change the sigma, this is the noise-free
    % version
    % for noise-case refer to Run01-OurApproaches_Noise
    sigma = 0;


    % the following coefficents are used in the arccosine-free version of the
    % algorithm
    a = fs/4/pi/sin(2*pi/N);
    b = 2*cos(2*pi/N);




    % %% filter
    % fnyq = fs/2;
    %
    %
    % w1 = f1/fnyq;
    % w2 = f2/fnyq;
    %
    % filterCoeffs = fir1(nFilter, [w1 w2]);

    %%






    t = 0:dt:tf;
    [row, col] = size(t);
    Nsamples = col;

    f = f0 *ones(size(t))+ sin(2*pi*1*t) + 0.5*sin(2*pi*6*t);
    %f = f0 *ones(size(t));
    %f = f0 *ones(size(t))+25*t-25*t.^2;
    f3 = 3*f;
    f5 = 5*f;

    f_max = max(f);
    f_min = min(f);

    w = 2*pi*f;

    w3 = 2*pi*f3;
    w5 = 2*pi*f5;




    phi0 = 0;

    theta1(1) = phi0;
    theta3(1) = phi0;
    theta5(1) = phi0;

    thetb1(1) = phi0+4*pi/3;
    thetb3(1) = phi0+4*pi/3;
    thetb5(1) = phi0+4*pi/3;

    thetc1(1) = phi0+2*pi/3;
    thetc3(1) = phi0+2*pi/3;
    thetc5(1) = phi0+2*pi/3;


    for k=2:length(t)
        theta1(k) = theta1(k-1)+w(k)*dt;
        theta3(k) = theta3(k-1)+w3(k)*dt;
        theta5(k) = theta5(k-1)+w5(k)*dt;

        thetb1(k) = thetb1(k-1)+w(k)*dt;
        thetb3(k) = thetb3(k-1)+w3(k)*dt;
        thetb5(k) = thetb5(k-1)+w5(k)*dt;

        thetc1(k) = thetc1(k-1)+w(k)*dt;
        thetc3(k) = thetc3(k-1)+w3(k)*dt;
        thetc5(k) = thetc5(k-1)+w5(k)*dt;
    end



    x = Am*cos(theta1);% + sigma*randn(size(t))+ h3Coeff*Am*cos(theta3)+ h5Coeff*Am*cos(theta5);
    y = Am*cos(thetb1);% + sigma*randn(size(t))+ h3Coeff*Am*cos(thetb3)+ h5Coeff*Am*cos(thetb5);
    z = Am*cos(thetc1);% + sigma*randn(size(t))+ h3Coeff*Am*cos(thetc3)+ h5Coeff*Am*cos(thetc5);

    [row, col] = size(x);
    power_DC = 0:col-1;

    %xh: x with harmonics
    %yh: y with harmonics
    %zh: z with harmonics
    xh = Am*cos(theta1) + A_DC_a*r_a.^power_DC + sigma*randn(size(t))+ h3Coeff*Am*cos(theta3)+ h5Coeff*Am*cos(theta5);
    yh = Am*cos(thetb1) + A_DC_b*r_b.^power_DC + sigma*randn(size(t))+ h3Coeff*Am*cos(thetb3)+ h5Coeff*Am*cos(thetb5);
    zh = Am*cos(thetc1) + A_DC_c*r_c.^power_DC + sigma*randn(size(t))+ h3Coeff*Am*cos(thetc3)+ h5Coeff*Am*cos(thetc5);

    % xhf: xh filtered
    % yhf: yh filtered
    % zhf: zh filtered

    % display('------------------------Removing DC using mimic filter (Benmouyal)--------------------')
    tau = 2*N;%2*T;%coeff_tau_a*T;
    temp1 = 1+tau-tau*cos(2*pi*f0/fs);
    temp2 = tau*sin(2*pi*f0/fs);
    den = temp1^2+temp2^2;
    den = sqrt(den);
    K_benmouyal = 1/den;
    filterCoeffs = K_benmouyal*[1+tau, -tau];
    %filterCoeffs = [1, -exp(-dt/tau)];
    nFilter = 2;
    xhf = filter(filterCoeffs,1,xh);
    yhf = filter(filterCoeffs,1,yh);
    zhf = filter(filterCoeffs,1,zh);


    % display('------------------------Removing DC using general case (My Approach)--------------------')
    % PS1_x = 0;
    % PS2_x = 0;
    % PS1_y = 0;
    % PS2_y = 0;
    % PS1_z = 0;
    % PS2_z = 0;
    % [row, col] = size(x);
    % for k=1:N
    %     PS1_x = PS1_x + xh(k);
    %     PS1_y = PS1_y + yh(k);
    %     PS1_z = PS1_z + zh(k);
    % end
    % PS2_x = PS1_x + xh(N+1) - xh(1);
    % PS2_y = PS1_y + yh(N+1) - yh(1);
    % PS2_z = PS1_z + zh(N+1) - zh(1);
    %
    % r_a_est = PS2_x/PS1_x;
    % r_b_est = PS2_y/PS1_y;
    % r_c_est = PS2_z/PS1_z;
    % A_DC_a_est = (r_a_est-1)/(r_a_est^N-1)*PS1_x;
    % A_DC_b_est = (r_b_est-1)/(r_b_est^N-1)*PS1_y;
    % A_DC_c_est = (r_c_est-1)/(r_c_est^N-1)*PS1_z;
    %
    % for k=1:N+1
    %     xhf(k) =  xh(k) - A_DC_a_est*r_a_est^(k-1);
    %     yhf(k) =  yh(k) - A_DC_b_est*r_b_est^(k-1);
    %     zhf(k) =  zh(k) - A_DC_c_est*r_c_est^(k-1);
    % end
    %
    % for k=N+2:col
    %     PS1_x = PS2_x;
    %     PS1_y = PS2_y;
    %     PS1_z = PS2_z;
    %
    %     PS2_x = PS1_x + xh(k) - xh(k-N);
    %     PS2_y = PS1_y + yh(k) - yh(k-N);
    %     PS2_z = PS1_z + zh(k) - zh(k-N);
    %
    %
    %     r_a_est = PS2_x/PS1_x;
    %     r_b_est = PS2_y/PS1_y;
    %     r_c_est = PS2_z/PS1_z;
    %     A_DC_a_est = (r_a_est-1)/(r_a_est^N-1)*PS1_x;
    %     A_DC_b_est = (r_b_est-1)/(r_b_est^N-1)*PS1_y;
    %     A_DC_c_est = (r_c_est-1)/(r_c_est^N-1)*PS1_z;
    %
    %     xhf(k) =  xh(k) - A_DC_a_est*r_a_est^N;
    %     yhf(k) =  yh(k) - A_DC_b_est*r_b_est^N;
    %     zhf(k) =  zh(k) - A_DC_c_est*r_c_est^N;
    %
    % end


    % display('------------------------Removing DC using general case (Recursive Removing)--------------------')
    % PS1_x = 0;
    % PS2_x = 0;
    % PS1_y = 0;
    % PS2_y = 0;
    % PS1_z = 0;
    % PS2_z = 0;
    % [row, col] = size(x);
    % for k=1:2:N-1
    %     PS1_x = PS1_x + xh(k);
    %     PS1_y = PS1_y + yh(k);
    %     PS1_z = PS1_z + zh(k);
    % end
    % for k=2:2:N
    %     PS2_x = PS2_x + xh(k);
    %     PS2_y = PS2_y + yh(k);
    %     PS2_z = PS2_z + zh(k);
    % end
    % r_a_est = PS2_x/PS1_x;
    % r_b_est = PS2_y/PS1_y;
    % r_c_est = PS2_z/PS1_z;
    % A_DC_a_est = (r_a_est^2-1)/r_a_est/(r_a_est^N-1)*PS1_x;
    % A_DC_b_est = (r_b_est^2-1)/r_b_est/(r_b_est^N-1)*PS1_y;
    % A_DC_c_est = (r_c_est^2-1)/r_c_est/(r_c_est^N-1)*PS1_z;
    %
    % for k=1:N
    %     xhf(k) =  xh(k) - A_DC_a_est*r_a_est^k;
    %     yhf(k) =  yh(k) - A_DC_b_est*r_b_est^k;
    %     zhf(k) =  zh(k) - A_DC_c_est*r_c_est^k;
    % end
    %
    % for k=N+1:col
    %     PS1_x_old = PS1_x;
    %     PS2_x_old = PS2_x;
    %     PS1_x = PS2_x_old;
    %     PS2_x = PS1_x_old + xh(k) - xh(k-N);
    %
    %
    %     PS1_y_old = PS1_y;
    %     PS2_y_old = PS2_y;
    %     PS1_y = PS2_y_old;
    %     PS2_y = PS1_y_old + yh(k) - yh(k-N);
    %
    %
    %     PS1_z_old = PS1_z;
    %     PS2_z_old = PS2_z;
    %     PS1_z = PS2_z_old;
    %     PS2_z = PS1_z_old + zh(k) - zh(k-N);
    %     r_a_est = PS2_x/PS1_x;
    %     r_b_est = PS2_y/PS1_y;
    %     r_c_est = PS2_z/PS1_z;
    %     A_DC_a_est = (r_a_est^2-1)/r_a_est/(r_a_est^N-1)*PS1_x;
    %     A_DC_b_est = (r_b_est^2-1)/r_b_est/(r_b_est^N-1)*PS1_y;
    %     A_DC_c_est = (r_c_est^2-1)/r_c_est/(r_c_est^N-1)*PS1_z;
    %
    %     xhf(k) =  xh(k) - A_DC_a_est*r_a_est^N;
    %     yhf(k) =  yh(k) - A_DC_b_est*r_b_est^N;
    %     zhf(k) =  zh(k) - A_DC_c_est*r_c_est^N;
    %
    % end

    % display('------------------------Removing DC using Simplified (Recursive Removing)--------------------')
    % PS1_x = 0;
    % PS2_x = 0;
    % PS1_y = 0;
    % PS2_y = 0;
    % PS1_z = 0;
    % PS2_z = 0;
    % [row, col] = size(x);
    % for k=1:2:N-1
    %     PS1_x = PS1_x + xh(k);
    %     PS1_y = PS1_y + yh(k);
    %     PS1_z = PS1_z + zh(k);
    % end
    % for k=2:2:N
    %     PS2_x = PS2_x + xh(k);
    %     PS2_y = PS2_y + yh(k);
    %     PS2_z = PS2_z + zh(k);
    % end
    % B1_x = 2/N/dt*(PS2_x-PS1_x);
    % B0_x = ((N+2)/N)*PS1_x-PS2_x;
    % B1_y = 2/N/dt*(PS2_y-PS1_y);
    % B0_y = ((N+2)/N)*PS1_y-PS2_y;
    % B1_z = 2/N/dt*(PS2_z-PS1_z);
    % B0_z = ((N+2)/N)*PS1_z-PS2_z;
    %
    % for k=1:N
    %     xhf(k) =  xh(k) - (B0_x+B1_x*k*dt);
    %     yhf(k) =  yh(k) - (B0_y+B1_y*k*dt);
    %     zhf(k) =  zh(k) - (B0_z+B1_z*k*dt);
    % end
    %
    % for k=N+1:col
    %     PS1_x_old = PS1_x;
    %     PS2_x_old = PS2_x;
    %     PS1_x = PS2_x_old;
    %     PS2_x = PS1_x_old + xh(k) - xh(k-N);
    %     B1_x = 2/N/dt*(PS2_x-PS1_x);
    %     B0_x = ((N+2)/N)*PS1_x-PS2_x;
    %
    %     PS1_y_old = PS1_y;
    %     PS2_y_old = PS2_y;
    %     PS1_y = PS2_y_old;
    %     PS2_y = PS1_y_old + yh(k) - yh(k-N);
    %     B1_y = 2/N/dt*(PS2_y-PS1_y);
    %     B0_y = ((N+2)/N)*PS1_y-PS2_y;
    %
    %     PS1_z_old = PS1_z;
    %     PS2_z_old = PS2_z;
    %     PS1_z = PS2_z_old;
    %     PS2_z = PS1_z_old + zh(k) - zh(k-N);
    %     B1_z = 2/N/dt*(PS2_z-PS1_z);
    %     B0_z = ((N+2)/N)*PS1_z-PS2_z;
    %
    %     xhf(k) =  xh(k) - (B0_x+B1_x*N*dt);
    %     yhf(k) =  yh(k) - (B0_y+B1_y*N*dt);
    %     zhf(k) =  zh(k) - (B0_z+B1_z*N*dt);
    % end


    % xhf = xh;%filter(filterCoeffs,1,xh);
    % yhf = yh;%filter(filterCoeffs,1,yh);
    % zhf = zh;%filter(filterCoeffs,1,zh);



    %% our approach for window length of 9

    M = 3;
    windowLength = 2*M+3; % = 9
    K = M; % K = (windowLength-3)/2

    firstIndex = K+2;

    for k=firstIndex:Nsamples-K-1
        Xk = x(k-K:k+K)';
        Xk_minus = x(k-K-1:k+K-1)';
        Xk_plus = x(k-K+1:k+K+1)';
        f_hatWin9(k) = fs/2/pi*acos(Xk'*(Xk_minus+Xk_plus)/2/(Xk'*Xk));
        f_hatWin9_arccosinefree(k) = f0 + a*(Xk'*(b*Xk-Xk_minus-Xk_plus))/(Xk'*Xk);

        Yk = y(k-K:k+K)';
        Yk_minus = y(k-K-1:k+K-1)';
        Yk_plus = y(k-K+1:k+K+1)';

        Zk = z(k-K:k+K)';
        Zk_minus = z(k-K-1:k+K-1)';
        Zk_plus = z(k-K+1:k+K+1)';


        numArg = Xk'*(Xk_minus+Xk_plus);
        numArg = numArg + Yk'*(Yk_minus+Yk_plus);
        numArg = numArg + Zk'*(Zk_minus+Zk_plus);

        denArg = 2*(Xk'*Xk);
        denArg = denArg + 2*(Yk'*Yk);
        denArg = denArg + 2*(Zk'*Zk);

        arg = numArg/denArg;
        f_hatWin9_3phase(k) = fs/2/pi*acos(arg);

        % Antonio Lopez Algorithm
        Xk = [x(k-3); x(k); x(k+3)];
        Xk_minus = [x(k-4); x(k-1); x(k+2)];
        Xk_plus = [x(k-2); x(k+1); x(k+4)];
        f_hat_lopez(k) = fs/2/pi*acos(Xk'*(Xk_minus+Xk_plus)/2/(Xk'*Xk));
    end

    [r, c] = size(t);
    Nsamples = c;

    Kmax =6;

    f_hat_1phase = f0*ones(Kmax+1,Kmax+1+nFilter);
    f_hat_3phase = f0*ones(Kmax+1,Kmax+1+nFilter);

    f_hat_1phase_withHarmonics = f0*ones(Kmax+1,Kmax+1+nFilter);
    f_hat_3phase_withHarmonics = f0*ones(Kmax+1,Kmax+1+nFilter);

    f_hat_1phase_withHarmonics_filtered = f0*ones(Kmax+1,Kmax+1+nFilter);
    f_hat_3phase_withHarmonics_filtered = f0*ones(Kmax+1,Kmax+1+nFilter);

    for K = 0:Kmax
        % K = 0,           Window Length = 3
        % K = 1,           Window Length = 5
        % K = 2,           Window Length = 7
        % K = 3,           Window Length = 9
        % K = 4,           Window Length = 11
        lengthWindow(K+1)  = 2*K+3;
        firstIndex = K+2+nFilter;% the first element of Xk_minus must be available: it must be the first element of the signal
        % k-K-1 = 1  ==>  k=K+2
        lastIndex = Nsamples-K-1;% the last element of Xk_plus must be available: it must be the first element of the signal
        % k+K+1 = Nsamples  ==> k = Nsamples-K-1
        for k=firstIndex:lastIndex
            Xk = x(k-K:k+K)';
            Xk_minus = x(k-K-1:k+K-1)';
            Xk_plus = x(k-K+1:k+K+1)';
            arg = Xk'*(Xk_minus+Xk_plus)/2/(Xk'*Xk);
            if(abs(arg)>1)
                f_hat_1phase(K+1,k) = f_hat_1phase(K+1,k-1);
            else
                f_hat_1phase(K+1,k) = fs/2/pi*acos(arg);
            end
            f_hat_1phase_arccosinefree(K+1,k) = f0 + a*(Xk'*(b*Xk-Xk_minus-Xk_plus))/(Xk'*Xk);

            Xkh = xh(k-K:k+K)';
            Xkh_minus = xh(k-K-1:k+K-1)';
            Xkh_plus = xh(k-K+1:k+K+1)';
            arg = Xkh'*(Xkh_minus+Xkh_plus)/2/(Xkh'*Xkh);
            if(abs(arg)>1)
                f_hat_1phase_withHarmonics(K+1,k) = f_hat_1phase_withHarmonics(K+1,k-1);
            else
                f_hat_1phase_withHarmonics(K+1,k) = fs/2/pi*acos(arg);
            end
            f_hat_1phase_arccosinefree_withHarmonics(K+1,k) = f0 + a*(Xkh'*(b*Xkh-Xkh_minus-Xkh_plus))/(Xkh'*Xkh);

            Xkhf = xhf(k-K:k+K)';
            Xkhf_minus = xhf(k-K-1:k+K-1)';
            Xkhf_plus = xhf(k-K+1:k+K+1)';
            arg = Xkhf'*(Xkhf_minus+Xkhf_plus)/2/(Xkhf'*Xkhf);
            if(abs(arg)>1)
                f_hat_1phase_withHarmonics_filtered(K+1,k) = f_hat_1phase_withHarmonics_filtered(K+1,k-1);
            else
                f_hat_1phase_withHarmonics_filtered(K+1,k) = fs/2/pi*acos(arg);
            end
            f_hat_1phase_arccosinefree_withHarmonics_filtered(K+1,k) = f0 + a*(Xkhf'*(b*Xkhf-Xkhf_minus-Xkhf_plus))/(Xkhf'*Xkhf);


            Yk = y(k-K:k+K)';
            Yk_minus = y(k-K-1:k+K-1)';
            Yk_plus = y(k-K+1:k+K+1)';

            Zk = z(k-K:k+K)';
            Zk_minus = z(k-K-1:k+K-1)';
            Zk_plus = z(k-K+1:k+K+1)';


            numArg = Xk'*(Xk_minus+Xk_plus);
            numArg = numArg + Yk'*(Yk_minus+Yk_plus);
            numArg = numArg + Zk'*(Zk_minus+Zk_plus);

            denArg = 2*(Xk'*Xk);
            denArg = denArg + 2*(Yk'*Yk);
            denArg = denArg + 2*(Zk'*Zk);

            arg = numArg/denArg;
            if(abs(arg)>1)
                f_hat_3phase(K+1,k) = f_hat_3phase(K+1,k-1);
            else
                f_hat_3phase(K+1,k) = fs/2/pi*acos(arg);
            end

            % three-phase arccosine-free
            num = a*(Xk'*(b*Xk-Xk_minus-Xk_plus));
            num = num + a*(Yk'*(b*Yk-Yk_minus-Yk_plus));
            num = num + a*(Zk'*(b*Zk-Zk_minus-Zk_plus));

            den = (Xk'*Xk);
            den = den + (Yk'*Yk);
            den = den + (Zk'*Zk);

            f_hat_3phase_arccosinefree(K+1,k) = f0 + num/den;

            % three-phase with harmonics
            Ykh = yh(k-K:k+K)';
            Ykh_minus = yh(k-K-1:k+K-1)';
            Ykh_plus = yh(k-K+1:k+K+1)';

            Zkh = zh(k-K:k+K)';
            Zkh_minus = zh(k-K-1:k+K-1)';
            Zkh_plus = zh(k-K+1:k+K+1)';


            numArg = Xkh'*(Xkh_minus+Xkh_plus);
            numArg = numArg + Ykh'*(Ykh_minus+Ykh_plus);
            numArg = numArg + Zkh'*(Zkh_minus+Zkh_plus);

            denArg = 2*(Xkh'*Xkh);
            denArg = denArg + 2*(Ykh'*Ykh);
            denArg = denArg + 2*(Zkh'*Zkh);

            arg = numArg/denArg;
            if(abs(arg)>1)
                f_hat_3phase_withHarmonics(K+1,k) = f_hat_3phase_withHarmonics(K+1,k-1);
            else
                f_hat_3phase_withHarmonics(K+1,k) = fs/2/pi*acos(arg);
            end

            % three-phase with harmonics arccosine-free
            num = a*(Xkh'*(b*Xkh-Xkh_minus-Xkh_plus));
            num = num + a*(Ykh'*(b*Ykh-Ykh_minus-Ykh_plus));
            num = num + a*(Zkh'*(b*Zkh-Zkh_minus-Zkh_plus));

            den = (Xkh'*Xkh);
            den = den + (Ykh'*Ykh);
            den = den + (Zkh'*Zkh);

            f_hat_3phase_arccosinefree_withHarmonics(K+1,k) = f0 + num/den;


            % three-phase with harmonics,  Filtered
            Ykhf = yhf(k-K:k+K)';
            Ykhf_minus = yhf(k-K-1:k+K-1)';
            Ykhf_plus = yhf(k-K+1:k+K+1)';

            Zkhf = zhf(k-K:k+K)';
            Zkhf_minus = zhf(k-K-1:k+K-1)';
            Zkhf_plus = zhf(k-K+1:k+K+1)';


            numArg = Xkhf'*(Xkhf_minus+Xkhf_plus);
            numArg = numArg + Ykhf'*(Ykhf_minus+Ykhf_plus);
            numArg = numArg + Zkhf'*(Zkhf_minus+Zkhf_plus);

            denArg = 2*(Xkhf'*Xkhf);
            denArg = denArg + 2*(Ykhf'*Ykhf);
            denArg = denArg + 2*(Zkhf'*Zkhf);

            arg = numArg/denArg;
            if(abs(arg)>1)
                f_hat_3phase_withHarmonics_filtered(K+1,k) = f_hat_3phase_withHarmonics_filtered(K+1,k-1);
            else
                f_hat_3phase_withHarmonics_filtered(K+1,k) = fs/2/pi*acos(arg);
            end

            % three-phase with harmonics arccosine-free, filtered
            num = a*(Xkhf'*(b*Xkhf-Xkhf_minus-Xkhf_plus));
            num = num + a*(Ykhf'*(b*Ykhf-Ykhf_minus-Ykhf_plus));
            num = num + a*(Zkhf'*(b*Zkhf-Zkhf_minus-Zkhf_plus));

            den = (Xkhf'*Xkhf);
            den = den + (Ykhf'*Ykhf);
            den = den + (Zkhf'*Zkhf);

            f_hat_3phase_arccosinefree_withHarmonics_filtered(K+1,k) = f0 + num/den;


        end
        %     figure()
        %     plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1)
        %     hold on
        %     plot(t(firstIndex:lastIndex), f_hat_1phase(K+1,firstIndex:lastIndex),'--','Color','blue','LineWidth',1)
        %     plot(t(firstIndex:lastIndex), f_hat_1phase_arccosinefree(K+1,firstIndex:lastIndex),'--','Color','green','LineWidth',1)
        %     plot(t(firstIndex:lastIndex), f_hat_3phase(K+1,firstIndex:lastIndex),'--','Color','red','LineWidth',1)
        %     plot(t(firstIndex:lastIndex), f_hat_3phase_arccosinefree(K+1,firstIndex:lastIndex),'--','Color','cyan','LineWidth',1)
        %     xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
        %     ylabel('f(t) and f_e_s_t(t))','FontSize', fontSize, 'FontWeight', 'bold')
        %     title(['data window length = ' int2str(2*K+3)],'FontSize', fontSize, 'FontWeight', 'bold')
        %     legend('exact frequency','estimated frequency')
        %     axis([0, tf, f_min-3, f_max+3])


        % with harmonics and/or noise
        se_1phase_withHarmonics = (f(firstIndex:lastIndex)-f_hat_1phase_withHarmonics(K+1,firstIndex:lastIndex)).^2;
        maxSE_1phase_withHarmonics = max(se_1phase_withHarmonics);
        MSE_1phase_withHarmonics(K+1) = mean(se_1phase_withHarmonics);
        MSE_1phase3phase_withHarmonics(K+1,1) = mean(se_1phase_withHarmonics);


        se_1phase_arccosinefree_withHarmonics = (f(firstIndex:lastIndex)-f_hat_1phase_arccosinefree_withHarmonics(K+1,firstIndex:lastIndex)).^2;
        maxSE_1phase_arccosinefree_withHarmonics = max(se_1phase_arccosinefree_withHarmonics);
        MSE_1phase_arccosinefree_withHarmonics(K+1) = mean(se_1phase_arccosinefree_withHarmonics);
        MSE_1phase3phase_withHarmonics(K+1,2) = mean(se_1phase_arccosinefree_withHarmonics);


        se_3phase_withHarmonics = (f(firstIndex:lastIndex)-f_hat_3phase_withHarmonics(K+1,firstIndex:lastIndex)).^2;
        maxSE_3phase_withHarmonics = max(se_3phase_withHarmonics);
        MSE_3phase_withHarmonics(K+1) = mean(se_3phase_withHarmonics);
        MSE_1phase3phase_withHarmonics(K+1,3) = mean(se_3phase_withHarmonics);


        se_3phase_arccosinefree_withHarmonics = (f(firstIndex:lastIndex)-f_hat_3phase_arccosinefree_withHarmonics(K+1,firstIndex:lastIndex)).^2;
        maxSE_3phase_arccosinefree_withHarmonics = max(se_3phase_arccosinefree_withHarmonics);
        MSE_3phase_arccosinefree_withHarmonics(K+1) = mean(se_3phase_arccosinefree_withHarmonics);
        MSE_1phase3phase_withHarmonics(K+1,4) = mean(se_3phase_arccosinefree_withHarmonics);

        % filtered
        se_1phase_withHarmonics_filtered = (f(firstIndex:lastIndex-nFilter/2)-f_hat_1phase_withHarmonics_filtered(K+1,firstIndex+nFilter/2:lastIndex)).^2;
        maxSE_1phase_withHarmonics_filtered = max(se_1phase_withHarmonics_filtered);
        MSE_1phase_withHarmonics_filtered(K+1) = mean(se_1phase_withHarmonics_filtered);
        MSE_1phase3phase_withHarmonics_filtered(K+1,1) = mean(se_1phase_withHarmonics_filtered);

        se_1phase_arccosinefree_withHarmonics_filtered = (f(firstIndex:lastIndex-nFilter/2)-f_hat_1phase_arccosinefree_withHarmonics_filtered(K+1,firstIndex+nFilter/2:lastIndex)).^2;
        maxSE_1phase_arccosinefree_withHarmonics_filtered = max(se_1phase_arccosinefree_withHarmonics_filtered);
        MSE_1phase_arccosinefree_withHarmonics_filtered(K+1) = mean(se_1phase_arccosinefree_withHarmonics_filtered);
        MSE_1phase3phase_withHarmonics_filtered(K+1,2) = mean(se_1phase_arccosinefree_withHarmonics_filtered);


        se_3phase_withHarmonics_filtered = (f(firstIndex:lastIndex-nFilter/2)-f_hat_3phase_withHarmonics_filtered(K+1,firstIndex+nFilter/2:lastIndex)).^2;
        maxSE_3phase_withHarmonics_filtered = max(se_3phase_withHarmonics_filtered);
        MSE_3phase_withHarmonics_filtered(K+1) = mean(se_3phase_withHarmonics_filtered);
        MSE_1phase3phase_withHarmonics_filtered(K+1,3) = mean(se_3phase_withHarmonics_filtered);

        se_3phase_arccosinefree_withHarmonics_filtered = (f(firstIndex:lastIndex-nFilter/2)-f_hat_3phase_arccosinefree_withHarmonics_filtered(K+1,firstIndex+nFilter/2:lastIndex)).^2;
        maxSE_3phase_arccosinefree_withHarmonics_filtered = max(se_3phase_arccosinefree_withHarmonics_filtered);
        MSE_3phase_arccosinefree_withHarmonics_filtered(K+1) = mean(se_3phase_arccosinefree_withHarmonics_filtered);
        MSE_1phase3phase_withHarmonics_filtered(K+1,4) = mean(se_3phase_arccosinefree_withHarmonics_filtered);


        %      figure()
        %      plot(t(firstIndex:lastIndex),mse,'Color','black','LineWidth',1)
        %      xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
        %      ylabel('MSE(f-f_e_s_t)','FontSize', fontSize, 'FontWeight', 'bold')
        %      title(['data window length = ' int2str(2*K+3)],'FontSize', fontSize, 'FontWeight', 'bold')
        %      axis([0, tf, 0, maxSE])
    end

    % for 1-phase case, the wibdow of length 3 is completely inefficent
    % therefore we ignore showing their MSE
    MSE_1phase3phase_withHarmonics(1,1)=0;
    MSE_1phase3phase_withHarmonics(1,2)=0;
    % figure(1)
    % bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
    % xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('mean(f-f_e_s_t)^2','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free', '3-phase', '3-phase arccosine-free')
    % title(['Comparison of four proposed approaches, 3^r^d harmonic = ' int2str(h3Coeff*100),'%, 5^t^h harmonic = ' int2str(h5Coeff*100),'%'],'FontSize', fontSize, 'FontWeight', 'bold')
    % %title(['Comparison of four proposed approaches, SNR = ' int2str(20*log10(1/sigma))],'FontSize', fontSize, 'FontWeight', 'bold')

    % for 1-phase case, the wibdow of length 3 is completely inefficent
    % therefore we ignore showing their MSE
    MSE_1phase3phase_withHarmonics_filtered(1,1)=0;
    MSE_1phase3phase_withHarmonics_filtered(1,2)=0;
    % figure(2)
    % bar(lengthWindow, MSE_1phase3phase_withHarmonics_filtered,'group')
    % xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('mean(f-f_e_s_t)^2','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free', '3-phase', '3-phase arccosine-free')
    % title(['Comparison of four proposed approaches, 3^r^d harmonic = ' int2str(h3Coeff*100),'%, 5^t^h harmonic = ' int2str(h5Coeff*100),'%, signal pre-filtered'],'FontSize', fontSize, 'FontWeight', 'bold')
    % %title(['Comparison of four proposed approaches, SNR = ' int2str(20*log10(1/sigma)) ', signal pre-filtered'],'FontSize', fontSize, 'FontWeight', 'bold')

    % figure(3)
    % plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1)
    % hold on
    % %plot(t(firstIndex:lastIndex), f_hat_1phase(4, firstIndex:lastIndex), 'red')
    % %plot(t(firstIndex:lastIndex), f_hat_1phase_withHarmonics_filtered(4,firstIndex:lastIndex), 'green')
    % plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), '--', 'Color','b','LineWidth',1)
    % plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), '--', 'Color','black','LineWidth',1)
    % plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), '--', 'Color','red','LineWidth',1)
    % plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), '--', 'Color','cyan','LineWidth',1)
    % xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
    % title('Comparison of exact and estimated frequency','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('exact frequency', '1Ph_WL9_DC', '1Ph_WL9_unDC', '3Ph_WL3_DC', '3Ph_WL3_unDC')

    % figure(4)
    % subplot(211)
    % bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
    % xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('mean(f-f_e_s_t)^2','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free', '3-phase', '3-phase arccosine-free')
    % title(['Comparison of four proposed approaches, 3^r^d harmonic = ' int2str(h3Coeff*100),'%, 5^t^h harmonic = ' int2str(h5Coeff*100),'%'],'FontSize', fontSize, 'FontWeight', 'bold')
    %
    % subplot(212)
    % bar(lengthWindow, MSE_1phase3phase_withHarmonics_filtered,'group')
    % xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('mean(f-f_e_s_t)^2','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free', '3-phase', '3-phase arccosine-free')
    % title(['Comparison of four proposed approaches, 3^r^d harmonic = ' int2str(h3Coeff*100),'%, 5^t^h harmonic = ' int2str(h5Coeff*100),'%, signal pre-filtered'],'FontSize', fontSize, 'FontWeight', 'bold')
    %
    [row, col] = size(MSE_1phase3phase_withHarmonics);

    for k=1:row
        MSE_1phase_New_withHarmonics(k, 1) = MSE_1phase3phase_withHarmonics(k, 1);
        MSE_1phase_New_withHarmonics(k, 2) = MSE_1phase3phase_withHarmonics(k, 2);
        MSE_3phase_New_withHarmonics(k, 1) = MSE_1phase3phase_withHarmonics(k, 3);
        MSE_3phase_New_withHarmonics(k, 2) = MSE_1phase3phase_withHarmonics(k, 4);

        MSE_1phase_New_withHarmonics_filtered(k, 1) = MSE_1phase3phase_withHarmonics_filtered(k, 1);
        MSE_1phase_New_withHarmonics_filtered(k, 2) = MSE_1phase3phase_withHarmonics_filtered(k, 2);
        MSE_3phase_New_withHarmonics_filtered(k, 1) = MSE_1phase3phase_withHarmonics_filtered(k, 3);
        MSE_3phase_New_withHarmonics_filtered(k, 2) = MSE_1phase3phase_withHarmonics_filtered(k, 4);
    end


    switch(iteration)
        case 1
            figure(1)
            subplot(411)
            bar(lengthWindow, MSE_1phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            legend('original version', 'arccosine-free version')
            title({'Decaying DC effect (only one phase polluted with decaying DC)';'1-phase case'},'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(412)
            bar(lengthWindow, MSE_3phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase', '3-phase arccosine-free')
            title(['3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')


            subplot(413)
            bar(lengthWindow, MSE_1phase_New_withHarmonics_filtered,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('1-phase (mimic filtered)', '1-phase arccosine-free (mimic filtered)')
            title(['1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(414)
            bar(lengthWindow, MSE_3phase_New_withHarmonics_filtered,'group')
            xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase (mimic filtered)', '3-phase arccosine-free (mimic filtered)')
            title(['3-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')

            figure(3)
            subplot(231)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 2, 'Color','k','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), '-.', 'Color','r','LineWidth',1)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), ':', 'MarkerSize', 3, 'Color','b','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), ':', 'Color','b','LineWidth',1)
%            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (only one phase polluted with decaying DC)','FontSize', fontSize, 'FontWeight', 'bold')
            title('Fig. 11a: only one phase polluted','FontSize', fontSize, 'FontWeight', 'bold')
%            legend('exact frequency', '1phase)')
            legend('exact', '1Ph')
%             legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 48, 54])

            subplot(232)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 2, 'Color','r','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), '-.', 'Color','r','LineWidth',1)
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 3, 'Color','k','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), ':', 'Color','b','LineWidth',1)
%            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
%            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (only one phase polluted with decaying DC)','FontSize', fontSize, 'FontWeight', 'bold')
            title('Fig. 11b: only one phase polluted','FontSize', fontSize, 'FontWeight', 'bold')
%            legend('exact frequency', '3phase')
            legend('exact', '3Ph')
%             legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 48, 54])


            subplot(234)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 2, 'Color','r','LineWidth',1.5)
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), '-.', 'Color','k','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), ':', 'MarkerSize', 3, 'Color','b','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), ':', 'Color','b','LineWidth',1.5)
            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (only one phase polluted with decaying DC)','FontSize', fontSize, 'FontWeight', 'bold')
            title({'Fig. 11d: only one phase polluted'; '(mimic filtered)'},'FontSize', fontSize, 'FontWeight', 'bold')
%            legend('exact frequency', '1phase+mimic)')
            legend('exact', '1Ph+mimic')
%             legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 50.17, 50.84])
            

            subplot(235)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 2, 'Color','r','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), '-.', 'Color','r','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), ':', 'MarkerSize', 3, 'Color','b','LineWidth',1.5)
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), '-.', 'Color','k','LineWidth',1.5)
            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
%            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (only one phase polluted with decaying DC)','FontSize', fontSize, 'FontWeight', 'bold')
            title({'Fig. 11e: only one phase polluted'; '(mimic filtered)'},'FontSize', fontSize, 'FontWeight', 'bold')
%            legend('exact frequency', '3phase+mimic')
            legend('exact', '3Ph+mimic')
%             legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 50.17, 50.84])

            
            figure(4)
            subplot(611)
            bar(lengthWindow, MSE_1phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            legend('original version', 'arccosine-free version')
            title({'Decaying DC effect';'Fig. 10a: 1-phase case'},'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(612)
            bar(lengthWindow, MSE_3phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase', '3-phase arccosine-free')
            title(['Fig. 10b: 3-phase case (only one phase polluted with decaying DC)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')


            subplot(614)
            bar(lengthWindow, MSE_1phase_New_withHarmonics_filtered,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('1-phase (mimic filtered)', '1-phase arccosine-free (mimic filtered)')
            title(['Fig. 10d: 1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(615)
            bar(lengthWindow, MSE_3phase_New_withHarmonics_filtered,'group')
%            xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase (mimic filtered)', '3-phase arccosine-free (mimic filtered)')
            title(['Fig. 10e: 3-phase case (mimic filtered) (only one phase polluted with decaying DC)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')



        case 2
            figure(2)
            subplot(411)
            bar(lengthWindow, MSE_1phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            legend('original version', 'arccosine-free version')
            title({'Decaying DC effect (all three phases polluted with the same decaying DC''s)';'1-phase case'},'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(412)
            bar(lengthWindow, MSE_3phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase', '3-phase arccosine-free')
            title(['3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')


            subplot(413)
            bar(lengthWindow, MSE_1phase_New_withHarmonics_filtered,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('1-phase (mimic filtered)', '1-phase arccosine-free (mimic filtered)')
            title(['1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 1-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')

            subplot(414)
            bar(lengthWindow, MSE_3phase_New_withHarmonics_filtered,'group')
            xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase (mimic filtered)', '3-phase arccosine-free (mimic filtered)')
            title('3-phase case (mimic filtered)','FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')

            figure(3)
            subplot(233)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), 'b-.', 'MarkerSize', 2, 'Color','r','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), 'b-.', 'Color','r','LineWidth',1)
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), '-.', 'MarkerSize', 3, 'Color','k','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), 'g:', 'Color','b','LineWidth',1)
%            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
%            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (all three phases polluted with the same decaying DC''s)','FontSize', fontSize, 'FontWeight', 'bold')
            title('Fig. 11c: three phases polluted identically','FontSize', fontSize, 'FontWeight', 'bold')
%            legend('exact frequency', '3phase')
            legend('exact', '3Ph')
%            legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 48, 54])

            subplot(236)
            plot(t(firstIndex:lastIndex), f(firstIndex:lastIndex),'Color','black','LineWidth',1.5)
            hold on
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics(4,firstIndex+nFilter/2:lastIndex), 'b-.', 'MarkerSize', 2, 'Color','r','LineWidth',1.5)
%            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_1phase_withHarmonics_filtered(4,firstIndex+nFilter/2:lastIndex), 'b-.', 'Color','r','LineWidth',1.5)
 %           plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics(1,firstIndex+nFilter/2:lastIndex), 'g:', 'MarkerSize', 3, 'Color','b','LineWidth',1.5)
            plot(t(firstIndex:lastIndex-nFilter/2), f_hat_3phase_withHarmonics_filtered(1,firstIndex+nFilter/2:lastIndex), '-.', 'Color','k','LineWidth',1.5)
            xlabel('time (sec)','FontSize', fontSize, 'FontWeight', 'bold')
%            ylabel('frequency (Hz)','FontSize', fontSize, 'FontWeight', 'bold')
%            title('Decaying DC effect on frequency estimation (all three phases polluted with the same decaying DC''s)','FontSize', fontSize, 'FontWeight', 'bold')
            title({'Fig. 11f: three phases polluted identically'; '(mimic filtered)'},'FontSize', fontSize, 'FontWeight', 'bold')
            legend('exact', '3Ph+mimic')
%            legend('exact frequency', '1phase (window length 9)', '1phase+mimic (window length 9)', '3phase (window length 3)', '3phase+mimic (window length 3)')
            axis([t(1), t(end), 50.17, 50.84])

            
            figure(4)

            subplot(613)
            bar(lengthWindow, MSE_3phase_New_withHarmonics,'group')
            %xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase', '3-phase arccosine-free')
            title(['Fig. 10c: 3-phase case (all three phases polluted with the same decaying DC''s)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case'],'FontSize', fontSize, 'FontWeight', 'bold')


            subplot(616)
            bar(lengthWindow, MSE_3phase_New_withHarmonics_filtered,'group')
            xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
            ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
            %legend('3-phase (mimic filtered)', '3-phase arccosine-free (mimic filtered)')
            title(['Fig. 10f: 3-phase case (mimic filtered) (all three phases polluted with the same decaying DC''s)'],'FontSize', fontSize, 'FontWeight', 'bold')
            %title(['Decaying DC effect on 3-phase case (mimic filtered)'],'FontSize', fontSize, 'FontWeight', 'bold')
            
    end

    % plot(samplingFreq, MSE_1PhaseWin9(1:15),'Color','black','LineWidth',1.5)
    % hold on
    % plot(samplingFreq , MSE_1PhaseArccosinefreeWin9(1:15), '-.','Color','black','LineWidth',1.5)
    % ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
    % title('for noisy signal, SNR = 60dB','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free')


    % figure(1)
    % subplot(311)
    % bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
    % xlabel('length of data window','FontSize', fontSize, 'FontWeight', 'bold')
    % ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
    % legend('1-phase', '1-phase arccosine-free', '3-phase', '3-phase arccosine-free')
    % title(['Fig. 9a: Comparison of four proposed approaches, with ' int2str(hInt1Coeff*100),'% of interharmonic of ', num2str(fInt1(1)), ' Hz'],'FontSize', fontSize, 'FontWeight', 'bold')
end