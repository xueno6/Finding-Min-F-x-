function result=equation(x)
T1 = 20*60;T2 = 50*60;T3 = 80*60;T4 = 140*60;NB1 = 20;NB2 = 15;NB3 = 15;NB4 = 12;
T = ((2+4*7)*3600 + T1*NB1 + T2*NB2 + T3*NB3 + T4*NB4)/952;   % //Overall time. The coming line giving the definition of T (cf # - General parameters) has to be commented and replaced by this one.
%// # - General parameters------------------------------------------------------ 
%// Definition of T (time) and R (radius)
%// T = 20;     // 1 T unit <=> 952 sec (15.873 min)
R = 1;     % // Since the equations are nondimensionalized with R as a reference space unit, our nondimensionalized bead' radius equals to one
%// "Size" of matrices - the real matrices have the size (J+1,M+1)
J = 5;
M = 2*floor(2*(T*J*J)/(R*R));    % //Based on the Courant–Friedrichs–Lewy condition

%// Deltat and Deltar - time & space steps for the numerical simulation
Deltat = T/M;
Deltar = R/J;

%// discrete instant t_m = (m-1) Deltat where Deltat = T / M
%// discrete position r_j = (j-1) Deltax where Deltar = R / J


%// Fitting parameters---------------------------------------------------------- /!\ TO BE MODIFIED WHEN IMPROVING THE FITTING /!\

%// Nondimensionalized diffusion coefficients taking water diffusion coefficient as a reference unit
CH2O =1;
CHO = x(1);%4e-1;
CHCO3 =x(2);% 1.2e-1;
CCO3 = x(3); % 1.5e-1;

%// K1, is the chemical equilibrium constant for the reaction: CO2 + HO- = HCO3-
%// K1 is supposed to be constant here.
K1 = x(4);%4e9;

%// ka is the Henry constant. It is the link between the partial pressure of water and the water
%// quantity on the surface: P(H2O) = ka * n(H2O)
ka = x(5);%1e-5;

%//kC is the ratio Volume(Bead)/Volume(Air) that will link P(CO2) with NHCO3 + NCO3
kC =x(6);% 4e-3;

%// K2 is the chemical equilibrium constant for the reaction: HCO3- + HO- = CO3-- + H2O. We'll make the assumption
%//that K2 is also a function of NH2O with K2 = kp * NH2O^p. This relation is true everywhere. Here we define the value
%// of kp and p.
kp =x(7);% 5e-12;
p = x(8);%9.5;

%// Peq is the reference pressure unit in our nondimensionlized equations
Peq = 5.104e7;% //unit = 0,1 Pa
%//Peq = n0.R.T in Pa with n0 = 1,94 mol/kg with 
%//a density d = 1,08 : 1m3 of bead <-> 1080 kg bead (# d' = 0,67 kg/L in the bucket) 
%//==> n0 = 2095 mol/m3 ==> Peq = n0.R.T = 5,104 10^6 Pa


 
 
 
%//--------------------------ii. Initial Conditions------------------------------

%//Initial condions outside the bead

PH2O = 15000 / Peq;%                      //initial PH2O: 15 000 / Peq means "150 000 Pa Nondimensionalized" since pressures are given in 0,1 Pa
%// In the case of oscillating PH2O (cf Boundary conditions below), this value will be used as a mean value or initial value for P(H2O).

PCO2 = 500 / Peq;%                        //initial PH2O: 500 / Peq means "5 000 Pa Nondimensionalized" since pressures are given in 0,1 Pa


%// We can calculate NH2O inside the bead (NH2Oins) thanks to Henry's law
NH2Oins = PH2O / ka;

%// Once you know nH2O, you can calculate K2 inside the beads (K2ins)
K2ins = kp * NH2Oins^p;
a4 = 0;
b4 = 1;

fa4 = 2 * ((ka * K2ins)/(PH2O * PCO2 * K1)) * a4*a4 + a4*(1 + 1/(PCO2 * K1)) - 1;
fb4 = 2 * ((ka * K2ins)/(PH2O * PCO2 * K1)) * b4*b4 + b4*(1 + 1/(PCO2 * K1)) - 1;

eps = 1e-9;
while (b4-a4) >eps
    X4 = (a4+b4)/2;
    fX4 = 2 * ((ka * K2ins)/(PH2O * PCO2 * K1)) * X4*X4 + X4*(1 + 1/(PCO2 * K1)) - 1;
    if (sign(fX4) == sign(fa4))
        a4 = X4;
        fa4 = fX4;
    else
        b4 = X4;
        fb4 = fX4;
    end
end

NHCO3ins = X4;

%// Once NHCO3 is known, it is simple to find NHO and NCO3.
NHOins = NHCO3ins / (K1 * PCO2);
NCO3ins = NHCO3ins * NHOins * ka * K2ins /PH2O;

%// Definition and initialization of the matrices
X = zeros(J,M);PCO2R = zeros(1,M+1);PH2OR = zeros(1,M+1);PCO2_OUT = zeros(M+1,1);PH2O_OUT = zeros(M+1,1);

MCO3 = zeros(J+1,M+1);MHCO3 = zeros(J+1,M+1);MHO = zeros(J+1,M+1);MH2O = zeros(J+1,M+1);

MCO3d = zeros(J+1,M);MHCO3d = zeros(J+1,M);MHOd = zeros(J+1,M);MH2Od = zeros(J+1,M);

MCO3dd = zeros(J+1,M);MHCO3dd = zeros(J+1,M);MHOdd = zeros(J+1,M);MH2Odd = zeros(J+1,M);

delta_electro = zeros(J+1,M);
MDIC = zeros(J+1,M+1);
K2 = zeros(J+1,M+1);
Vd = zeros(J+1,M);
Vdd = zeros(J+1,M);

JH2O = zeros(J+1,M);JHO = zeros(J+1,M);JCO3 = zeros(J+1,M);JHCO3 = zeros(J+1,M);JDIC = zeros(J+1,M);


%// Initial overall equilibrium everywhere in the bead. At the begining, the system is at the corresponding steady state
for j = 1:J+1
    MHO(j,1) = NHOins;
    MHCO3(j,1) = NHCO3ins;
    MCO3(j,1) = NCO3ins;
    MH2O(j,1) = NH2Oins;
    MDIC(j,1) = MCO3(j,1) + MHCO3(j,1);
end

for m = 1 : M+1
    PH2OR(m) = PH2O;
end 
m00 = round((2*3600*M)/(952*T))+1;
m0 = round((7*3600*M)/(952*T));
%//
%//
%////T1 = 20 min - 20 periods
w1 = 2*pi/(T1/952);
m1 = m00;
m11 = m1 + round((T1*NB1*M)/(952*T));
if NB1 > 0
for m = m1 : m11
    PH2OR(m) = PH2O + PH2O/3*sin(w1*(m-1-(m1-1))*T/M);
end 
end
% //
% //
% ////T2 = 50 min - 20 periods
w2 = 2*pi/(T2/952);
m2 = m11 + m0;
m22 = m2 + round((T2*NB2*(M))/(952*T));
if NB2 > 0 
for m = m2 : m22
    PH2OR(m) = PH2O + PH2O/3*sin(w2*(m-1-(m2-1))*T/M);
end 
end
% //
% ////T3 = 80 min - 15 periods
w3 = 2*pi/(T3/952);
m3 = m22 + m0;
m33 = m3 + round((T3*NB3*(M))/(952*T));
if NB3 > 0 
for m = m3 : m33
    PH2OR(m) = PH2O + PH2O/3*sin(w3*(m-1-(m3-1))*T/M);
end
end
% //
% //
% ////T4 = 140 min - 10 periods
w4 = 2*pi/(T4/952);
m4 = m33 + m0;
m44 = m4 + round((T4*NB4*(M))/(952*T));
if NB4 > 0 
for m = m4 : m44
    PH2OR(m) = PH2O + PH2O/3*sin(w4*(m-1-(m4-1))*T/M);
end 
end
%%------------------------------------------------------------
for m = 1:M+1
    %// By definition for H2O
    NH2OR(m) = PH2OR(m) / ka;
    %// Once you know NH2O, you can calculate K2
    K2R(m) = kp * NH2OR(m).^p;
end
%// We also know PCO2 at the surface (PCO2R) at the begining for m = 1 (t_1)
PCO2R(1) = PCO2;

for m = 1:M
    for j = 2:J
        %// Discrete first derivatives
        MHOd(j,m)=(MHO(j+1,m) - MHO(j-1,m))/(2*Deltar);
        MHCO3d(j,m)=(MHCO3(j+1,m) - MHCO3(j-1,m))/(2*Deltar);
        MCO3d(j,m)=(MCO3(j+1,m) - MCO3(j-1,m))/(2*Deltar);
        MH2Od(j,m)=(MH2O(j+1,m) - MH2O(j-1,m))/(2*Deltar);
        Vd(j,m)=(CHO*MHOd(j,m)+CHCO3*MHCO3d(j,m)+2*CCO3*MCO3d(j,m))/(CHO*MHO(j,m)+CHCO3*MHCO3(j,m)+4*CCO3*MCO3(j,m));
        
        %// Discrete second derivatives
        MHOdd(j,m)=(MHO(j+1,m) - 2*MHO(j,m) + MHO(j-1,m))/(Deltar * Deltar);
        MHCO3dd(j,m)=(MHCO3(j+1,m) - 2*MHCO3(j,m) + MHCO3(j-1,m))/(Deltar * Deltar);
        MCO3dd(j,m)=(MCO3(j+1,m) - 2*MCO3(j,m) + MCO3(j-1,m))/(Deltar * Deltar);
        MH2Odd(j,m)=(MH2O(j+1,m) - 2*MH2O(j,m) + MH2O(j-1,m))/(Deltar * Deltar);
        Vdd(j,m)=(CHO*MHOdd(j,m)+CHCO3*MHCO3dd(j,m)+2*CCO3*MCO3dd(j,m))/(CHO*MHO(j,m)+CHCO3*MHCO3(j,m)+4*CCO3*MCO3(j,m))-(CHO*MHOd(j,m)+CHCO3*MHCO3d(j,m)+2*CCO3*MCO3d(j,m))*(CHO*MHOd(j,m)+CHCO3*MHCO3d(j,m)+4*CCO3*MCO3d(j,m))/((CHO*MHO(j,m)+CHCO3*MHCO3(j,m)+4*CCO3*MCO3(j,m))^2);

        DeltaNHO2 = Deltat*CHO/9*(MHOdd(j,m) - MHOd(j,m)*Vd(j,m) - MHO(j,m)*Vdd(j,m) + 2/((j-1)*Deltar)*MHOd(j,m) - 2/((j-1)*Deltar)*MHO(j,m)*Vd(j,m));
        DeltaNHCO32 = Deltat*CHCO3/9*(MHCO3dd(j,m) - MHCO3d(j,m)*Vd(j,m) - MHCO3(j,m)*Vdd(j,m) + 2/((j-1)*Deltar)*MHCO3d(j,m) - 4/((j-1)*Deltar)*MHCO3(j,m)*Vd(j,m));
        DeltaNCO32 = Deltat*CCO3/9*(MCO3dd(j,m) - 2*MCO3d(j,m)*Vd(j,m) - 2*MCO3(j,m)*Vdd(j,m) + 2/((j-1)*Deltar)*MCO3d(j,m) - 2/((j-1)*Deltar)*MCO3(j,m)*Vd(j,m));
        DeltaNH2O2 = Deltat/9*(MH2Odd(j,m) + 2/((j-1)*Deltar)*MH2Od(j,m));
        
%         // Alpha (X)
        if (MHO(j,m) + DeltaNHO2) > (MHCO3(j,m) + DeltaNHCO32)
            a = -(MHCO3(j,m) + DeltaNHCO32);
        else
            a = - (MHO(j,m) + DeltaNHO2);
        end

        if (MCO3(j,m) + DeltaNCO32) > (MH2O(j,m) + DeltaNH2O2)
            b = (MH2O(j,m) + DeltaNH2O2);
        else
            b = (MCO3(j,m) + DeltaNCO32);
        end

        %// X must be comprised between a and b : a<X<b
        a2 = a;
        b2 = b;

        K2a = kp * (MH2O(j,m) + DeltaNH2O2 - a)^p;
        K2b = kp * (MH2O(j,m) + DeltaNH2O2 - b)^p;

        fa = K2a * (MHO(j,m) + DeltaNHO2 + a) * (MHCO3(j,m) + DeltaNHCO32 + a) - (MCO3(j,m) + DeltaNCO32 - a) * (MH2O(j,m) + DeltaNH2O2 - a);
        fb = K2b * (MHO(j,m) + DeltaNHO2 + b) * (MHCO3(j,m) + DeltaNHCO32 + b) - (MCO3(j,m) + DeltaNCO32 - b) * (MH2O(j,m) + DeltaNH2O2 - b);

        if (b2-a2)<0 
            warning ("absurd alpha boundaries")
            result=-1;
            return
        end

        eps = 1e-9;
        while (b2-a2) >eps
            if (fa*fb>0) 
                warning("fa and fb have the same sign")
                result=-1;
                return
            end
            X(j,m) = (a2+b2)/2;
            K2X = kp * (MH2O(j,m) + DeltaNH2O2 - X(j,m))^p;
            fX = K2X * (MHO(j,m) + DeltaNHO2 + X(j,m)) * (MHCO3(j,m) + DeltaNHCO32 + X(j,m)) - (MCO3(j,m) + DeltaNCO32 - X(j,m)) * (MH2O(j,m) + DeltaNH2O2 - X(j,m));
            if (sign(fX) == sign(fa)) 
                a2 = X(j,m);
                fa = fX;
                K2a = K2X;
            else
                b2 = X(j,m);
                fb = fX;
                K2b = K2X;
            end
        end

        if X(j,m) - a < 0 
            warning ("alpha is not bounded IN")
            result=-1;
            return
        end

        if X(j,m) - b > 0 
            warning ("alpha is not bounded IN")
            result=-1;
            return
        end


        %// At t + Deltat: calculation of the new concentration and flux at t_(m+1)
        MHO(j,m+1) = MHO(j,m) + DeltaNHO2 + X(j,m);
        MHCO3(j,m+1) = MHCO3(j,m) + DeltaNHCO32 + X(j,m);
        MCO3(j,m+1) = MCO3(j,m) + DeltaNCO32 - X(j,m);
        MH2O(j,m+1) = MH2O(j,m) + DeltaNH2O2 - X(j,m);
        MDIC(j,m+1) = MCO3(j,m+1) + MHCO3(j,m+1);

        JHO(j,m)= - CHO * MHOd(j,m) + CHO * Vd(j,m) * MHO(j,m);
        JH2O(j,m)= - CH2O * MH2Od(j,m);
        JHCO3(j,m)= - CHCO3 * MHCO3d(j,m) + CHCO3 * Vd(j,m) * MHCO3(j,m);
        JCO3(j,m)= - CCO3 * MCO3d(j,m) + 2 * CCO3 * Vd(j,m) * MCO3(j,m);
        JDIC(j,m)= JHCO3(j,m) + JCO3(j,m);
    end


    MHOd(1,m)=0;
    MHCO3d(1,m)=0;
    MCO3d(1,m)=0;
    MH2Od(1,m)=0;
    Vd(1,m)=0;

        %// Discrete second derivatives
    MHOdd(1,m)=2*(MHO(2,m) - MHO(1,m))/(Deltar * Deltar);
    MHCO3dd(1,m)=2*(MHCO3(2,m) - MHCO3(1,m))/(Deltar * Deltar);
    MCO3dd(1,m)=2*(MCO3(2,m) - MCO3(1,m))/(Deltar * Deltar);
    MH2Odd(1,m)=2*(MH2O(2,m) - MH2O(1,m))/(Deltar * Deltar);
    Vdd(1,m)=(CHO*MHOdd(1,m)+CHCO3*MHCO3dd(1,m)+2*CCO3*MCO3dd(1,m))/(CHO*MHO(1,m)+CHCO3*MHCO3(1,m)+4*CCO3*MCO3(1,m));  


%         // Quantities diffusing between t_m and t_(m+1) at the center
%         // With the charged species, we consider that the flux is the sum of the
%         // diffusion current and the drift current.
        DeltaNHO2 = Deltat*CHO/9*(MHOdd(1,m) - MHO(1,m)*Vdd(1,m));
        DeltaNHCO32 = Deltat*CHCO3/9*(MHCO3dd(1,m) - MHCO3(1,m)*Vdd(1,m));
        DeltaNCO32 = Deltat*CCO3/9*(MCO3dd(1,m) - 2*MCO3(1,m)*Vdd(1,m));
        DeltaNH2O2 = Deltat/9*MH2Odd(1,m);


    %// Alpha (X)
    if (MHO(1,m) + DeltaNHO2) > (MHCO3(1,m) + DeltaNHCO32) 
        a = -(MHCO3(1,m) + DeltaNHCO32);
    else
        a = - (MHO(1,m) + DeltaNHO2);
    end

    if (MCO3(1,m) + DeltaNCO32) > (MH2O(1,m) + DeltaNH2O2) 
        b = (MH2O(1,m) + DeltaNH2O2);
    else
        b = (MCO3(1,m) + DeltaNCO32);
    end

    %// X must be comprised between a and b : a<X<b
    a2 = a;
    b2 = b;

    K2a = kp * (MH2O(1,m) + DeltaNH2O2 - a)^p;
    K2b = kp * (MH2O(1,m) + DeltaNH2O2 - b)^p;

    fa = K2a * (MHO(1,m) + DeltaNHO2 + a) * (MHCO3(1,m) + DeltaNHCO32 + a) - (MCO3(1,m) + DeltaNCO32 - a) * (MH2O(1,m) + DeltaNH2O2 - a);
    fb = K2b * (MHO(1,m) + DeltaNHO2 + b) * (MHCO3(1,m) + DeltaNHCO32 + b) - (MCO3(1,m) + DeltaNCO32 - b) * (MH2O(1,m) + DeltaNH2O2 - b);


    if (b2-a2)<0 
        warning ("absurd alpha boundaries")
        result=-1;
        return
    end

    eps = 1e-9;
    while (b2-a2) >eps
        if (fa*fb>0) 
            warning("fa and fb have the same sign")
            result=-1;
            return
        end
        X(1,m) = (a2+b2)/2;
        K2X = kp * (MH2O(1,m) + DeltaNH2O2 - X(1,m))^p;
        fX = K2X * (MHO(1,m) + DeltaNHO2 + X(1,m)) * (MHCO3(1,m) + DeltaNHCO32 + X(1,m)) - (MCO3(1,m) + DeltaNCO32 - X(1,m)) * (MH2O(1,m) + DeltaNH2O2 - X(1,m));
        if (sign(fX) == sign(fa)) 
            a2 = X(1,m);
            fa = fX;
            K2a = K2X;
        else
            b2 = X(1,m);
            fb = fX;
            K2b = K2X;
        end
    end

    if X(1,m) - a < 0 
        warning ("alpha is not bounded CENTER")
        result=-1;
        return
    end

    if X(1,m) - b > 0 
        warning ("alpha is not bounded CENTER")
        result=-1;
        return
    end
    
    
        %// At t + Deltat: calculation of the new concentration and flux at t_(m+1)
    MHO(1,m+1) = MHO(1,m) + DeltaNHO2 + X(1,m);
    MHCO3(1,m+1) = MHCO3(1,m) + DeltaNHCO32 + X(1,m);
    MCO3(1,m+1) = MCO3(1,m) + DeltaNCO32 - X(1,m);
    MH2O(1,m+1) = MH2O(1,m) + DeltaNH2O2 - X(1,m);
    MDIC(1,m+1) = MCO3(1,m+1) + MHCO3(1,m+1);

    JHO(1,m)= 0;
    JH2O(1,m)= 0;
    JHCO3(1,m)= 0;
    JCO3(1,m)= 0;
    JDIC(1,m)= 0;

     
                        if m > 1 
                            Delta_carbon = 0;
                            for j = 1 : J
                                Delta_carbon = Delta_carbon + ((j-0.5)*Deltar)^2*0.5*((MHCO3(j+1,m) + MHCO3(j,m)) - (MHCO3(j,m-1) + MHCO3(j+1,m-1)) + (MCO3(j+1,m) + MCO3(j,m)) - (MCO3(j+1,m-1) + MCO3(j,m-1)));%// we can approximate the bead to homogeneous nested shells with species' concentrations changing from one shell to the other.
                            end
                            PCO2R(m+1) =  PCO2R(m) - 3 * kC * Deltar * Delta_carbon;
                        else
                            PCO2R(m+1) =  PCO2R(m);
                        end 
                        JDIC(J+1,m) = 3*(PCO2R(m+1)-PCO2R(m))/(kC*Deltat); %// carbon flux at the boundary
                
                        if PCO2R(m+1)<0 
                            warning ("Negative CO2 partial Pressure")
                            result=-1;
                            return 
                        end
    
        a4 = 0;
        b4 = 1;
    
        fa4 = 2 * ((ka * K2R(m+1))/(PH2OR(m+1) * PCO2R(m+1) * K1)) * a4*a4 + a4*(1 + 1/(PCO2R(m+1) * K1)) - 1;
        fb4 = 2 * ((ka * K2R(m+1))/(PH2OR(m+1) * PCO2R(m+1) * K1)) * b4*b4 + b4*(1 + 1/(PCO2R(m+1) * K1)) - 1;
    
        eps = 1e-9;
        while (b4-a4) >eps
            X4 = (a4+b4)/2;
            fX4 = 2 * ((ka * K2R(m+1))/(PH2OR(m+1) * PCO2R(m+1) * K1)) * X4*X4 + X4*(1 + 1/(PCO2R(m+1) * K1)) - 1;
            if (sign(fX4) == sign(fa4)) 
                a4 = X4;
                fa4 = fX4;
            else
                b4 = X4;
                fb4 = fX4;
            end
        end
    
        NHCO3R(m+1) = X4;
    
        %// Once NHCO3R is known, it is simple to find NHOR and NCO3R.
        NHOR(m+1) = NHCO3R(m+1) / (K1 * PCO2R(m+1));
        NCO3R(m+1) = NHCO3R(m+1) * NHOR(m+1) * ka * K2R(m+1) /PH2OR(m+1);


        %// New boundary conditions at r = R with r(J+1) = R
        MHCO3(J+1,m+1) = NHCO3R(m+1);
        MHO(J+1,m+1) = NHOR(m+1);
        MCO3(J+1,m+1) = NCO3R(m+1);
        MH2O(J+1,m+1) = NH2OR(m+1);
        MDIC(J+1,m+1) = MCO3(J+1,m+1) + MHCO3(J+1,m+1);

        %//calculating the flux at the borders: r(J+1) = R by doing a linear extrapolation of concentrations' first derivative
        MHOd(J+1,m)= MHOd(J,m) + MHOdd(J,m) * Deltar;
        MHCO3d(J+1,m)= MHCO3d(J,m) + MHCO3dd(J,m) * Deltar;
        MCO3d(J+1,m)= MCO3d(J,m) + MCO3dd(J,m) * Deltar;
        MH2Od(J+1,m)= MH2Od(J,m) + MH2Odd(J,m) * Deltar;
        Vd(J+1,m)= Vd(J,m) + Vdd(J,m) * Deltar;
        
        
        JHO(J+1,m)= - CHO * MHOd(J+1,m) + CHO * Vd(J+1,m) * MHO(J+1,m);
        JH2O(J+1,m)= - CH2O * MH2Od(J+1,m);
        JHCO3(J+1,m)= - CHCO3 * MHCO3d(J+1,m) + CHCO3 * Vd(J+1,m) * MHCO3(J+1,m);
        JCO3(J+1,m)= - CCO3 * MCO3d(J+1,m) + 2 * CCO3 * Vd(J+1,m) * MCO3(J+1,m);
end

delta_eletromax = 0;
for m=1:M
    for j = 1:J+1
        delta_electro(j,m) = (MHCO3(j,m+1) + MHO(j,m+1) + 2*MCO3(j,m+1)) - (MHO(j,m) + MHCO3(j,m) + 2*MCO3(j,m));
        if abs(delta_electro(j,m)) > delta_eletromax 
            delta_electromax = abs(delta_electro(j,m));
%             jmax = j;
%             mmax = m;
        end
    end
end
% //--------------------------------------------- Carbon conservation
% // We want to check the carbon conservation during the whole process and the evolution of all the species

for m=1:M+1
    H2Otot(m) = 0;
    for j = 1:J
        H2Otot(m) = H2Otot(m) + ((j-0.5)*Deltar)^2*0.5*(MH2O(j+1,m)+MH2O(j,m));
    end
end


for m=1:M+1
    HOtot(m) = 0;
    for j = 1:J
        HOtot(m) = HOtot(m) + ((j-0.5)*Deltar)^2*0.5*(MHO(j+1,m)+MHO(j,m));
    end
end


for m=1:M+1
    HCO3tot(m) = 0;
    for j = 1:J
        HCO3tot(m) = HCO3tot(m) + ((j-0.5)*Deltar)^2*0.5*(MHCO3(j+1,m)+MHCO3(j,m));
    end
end

for m=1:M+1
    CO3tot(m) = 0;
    for j = 1:J
        CO3tot(m) = CO3tot(m) + ((j-0.5)*Deltar)^2*0.5*(MCO3(j+1,m)+MCO3(j,m));
    end
end

for m =1 : M+1
carbon_quantity(m) = 0;
end

absolute_carbon_var = 0;
if kC > 0 
    for m=1:M+1
        carbon_quantity(m) = PCO2R(m) + 3*kC*Deltar*(HCO3tot(m) + CO3tot(m));
        if m >1 
            absolute_carbon_var = absolute_carbon_var + abs(carbon_quantity(m)-carbon_quantity(m-1));
        end
    end
end

for m =1 : M+1
PCO2_OUT(m,1) = PCO2R(1,m)*Peq;
PH2O_OUT(m,1) = PH2OR(1,m)*Peq;
end

% //------------------------------vii. Plots -------------------------------------
input = (PH2OR'*Peq)';
output = (PCO2R'*Peq)';
result = calculateDistance([input; output], 24000, 0, @evaluateC);
end





