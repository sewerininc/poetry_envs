% read_simulation.m %%%%%%%%%
%%           Np     Ee  beta    Eion    Mion    Drift   F3u     F2u     F1u     EXTu    saddle  Center  EXTd    U_Fd1    U_Fd2  U_Fd3    U_M1    U_M2         % %
%%%                                                                                                                                                     % %
%%% TQual = 3
% Elect_1     200   2.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %        
% Elect_2     2000  2.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %        
% Elect_3     2000  1.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %        
% Elect_4     2000  3.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_5     2000  3.5 2.0     3000    50      25      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_6     2000  3.5 2.0     3000    50      35      0       350     0       500     0       50      500     0       350     0        -100     0        % %
%%% TQual = 150 to get better resloution on the dissociation point, but takes longer to run %
% Elect_7     2000  3.5 2.0     3000    50      35      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_8     2000  3.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_9     2000  1.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_10    2000  0.5 2.0     3000    50      50      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_11    2000  0.5 2.0     3000    50      75      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_12    2000  1.5 2.0     3000    50      75      0       350     0       500     0       50      500     0       350     0        -100     0        % %
% Elect_13    2000  1.5 2.0     3000    50     100      0       350     0       500     0       50      500     0       350     0        -100     0        % %

%Simulations for minimising the drift potential, and looking for potential
%setups for measures with less off-set for the energy of the
%photoelectrons. 
% Mikkel_02   2000  1.5 2.0     3000    50     5        0        5      0       20      0       50      20      0       10      0         10     10     By=50 gauss   % %
% Mikkel_04   2000  1.5 2.0     3000    50     3        0        5      0       15      0       50      15      0       10      0         10     10     By=50 gauss   % %
% Mikkel_05   2000  1.5 2.0     3000    50     1        0        5      0       15      0       50      15      0       10      0         10     10     By=50 gauss   % %
% Mikkel_06   2000  1.5 2.0     3000    50     2        0        5      0       10      0       50      10      0       10      0         10     10     By=100 gauss   % %
%Photoelectrons hits the detector well for the upgoing trajectory, but
%there is minor disturbences for the downwards path for 06. 
% Mikkel_07   2000  1.5 2.0     3000    50     1        0        15      0       5      0       50      5      0       15      0         15     15     By=100 gauss   % %
%this shows low disturbence in before the B-field is turned on, the low Ext
%is compensated by a slightly higher f2

%energy is 0.861, resembling oxygen, and beta is -1. 
% Mikkel_8    2000  1.5 2.0     3000    50     5        0        20      0       10      0       50      10      0       20      0         20     20     By=100 gauss   % %
% Mikkel_9    2000  1.5 2.0     3000    50     5        0        20      0       10      0       50      10      0       20      0         20     20     By=100 gauss   % %
% Mikkel_10    2000  1.5 2.0     3000    50     5        0        20      0       10      0       50      10      0       20      0         20     20     By=0 gauss   % %
clear

%files = [  3   4   5   6  7    8  9  10 11]
files = [243]
%files = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68]
pot   = [-10 -20   0  -5 -2  -1 -3  -4 -6]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program='sim8';

do_load_simulation_files=1;           % 0=NO, 1-Yes
    file_directory= 'O:\Nat_IFA-xring\CROSSED_BEAMS\SIMULATION_ELECTRON\'
    file_name='MMRUN_'          % full filename is: file_name#.dat
        n_start = 243                 % start at file no. 262
        n_stop  = 243;                 % stop at file no. 266
    %file_name='photon_2021_'          % full filename is: file_name#.dat
    %    n_start = 262                 % start at file no. 262
    %    n_stop  = 266;                 % stop at file no. 266

V_bias=00

do_define_special_inputs=1;         %     
        Xmid=75*1+0.0;              % mid point in X-axis (mm)
        Ymid=150*0;                 % mid point in Y-axis (mm)
        Zmid=75*0;                  % mid point in Z-axis (mm)
        Tmid=0.0;
 do_assignment_of_data=1            %     
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Natural constants  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c Qe AMU me ao to kBoltz au2eV eV2au
        c= 2.99792458e8             %2.99792558e8;  % m/s
        Qe=1.602176634e-19          %1.6022e-19;    % e
        AMU=931.4908918e6/(c^2)*Qe  % eV s2/m2
        me=1822.8885;               % AMU in units of the electron mass
        ao=0.529177;                % atomic length unit (Angstroem)
        to=2.41889e-17;             % atomic time unit (s) 
        kBoltz=3.168140552e-6;      % Boltzmann constant 
        au2eV=27.2113961;
        eV2au=1/au2eV;
        eV2cm=8065.5409;
        cm2eV=1/eV2cm;
        mm_microsec2m_sec=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L_DET1 = 0.263;
L_DET2 = 0.967;
L_DET3 = 0.975;

do_HEX=0;
do_DLD=0;
do_HCH=1;
do_ELEC=0;
% %%%%
% Vbias=500
% M=26*AMU;
% E=(4200-Vbias)*Qe
% v=sqrt(2*E/M)
% TOF_DLD=L_DET2/v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(do_load_simulation_files==1)
%%% Loading simulation files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear NumberOfRecordedParamters NumberOfInfoLines SimulationInfo
    NumberOfRecordedParamters=10;
    NumberOfInfoLines=3;
    SimulationInfo=zeros(NumberOfInfoLines,NumberOfRecordedParamters);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear FILE_DATA ALL_DATA ic
        ALL_DATA=zeros(1,NumberOfRecordedParamters);
        ic=0;
        for n_file=n_start:n_stop
           'loading file number'
            %n_file=105
            data_file=strcat(file_directory,file_name,num2str(n_file),'.dat');
                fid=fopen(data_file,'rt');
                switch program
                    case ('sim7')
                    [SimulationInfo,count]=fscanf(fid,'%20f',[NumberOfRecordedParamters,NumberOfInfoLines]);
                    fgets(fid);
                    Textline=fgets(fid);
                    case ('sim8')
                        
                    
                    Textline=fgets(fid);
                    %fgets(fid);
                    [SimulationInfo,count]=fscanf(fid,'%20f',[NumberOfRecordedParamters,NumberOfInfoLines]);
                end
        
                [FILE_DATA,count]=fscanf(fid,'%20f',[NumberOfRecordedParamters,inf]);
                fclose(fid);
             
                 a=size(FILE_DATA);
                a_2 = a(2)
                file_no = zeros(1,a_2)+n_file;
                FILE_DATA = [FILE_DATA; file_no];
                
             if (ic==0)
                ALL_DATA=FILE_DATA';
                ic=1;
             else    
                ALL_DATA=[ALL_DATA; FILE_DATA'];
             end
      end % Loading simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %do_load_simulation_files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(do_assignment_of_data==1)
%%%% Asignment of the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear n_ion_index TOF_index Mass_index Charge_index
clear X_index Y_index Z_index
clear Vx_index Vy_index Vz_index
clear row_all col_all

[row_all,col_all]=size(ALL_DATA)
    number_index=1;
    TOF_index=2;
    Mass_index=3;
    Charge_index=4;
    X_index=5;
    Y_index=6;
    Z_index=7;
    Vx_index=8;
    Vy_index=9;
    Vz_index=10;
    file_index =11;
 clear nr   
 nr=8;             % depends on setting nr=2 times the number of lines for each ion
 %%% Ions starting parameters %%%%%%%%%%%%%%%%%%%%%%%%
 clear Ion_n Ion_TOF Ion_Mass Ion_Charge 
 clear Ion_X  Ion_Y Ion_Z 
 clear Ion_Vx Ion_Vy Ion_Vz
 Ion_n      = ALL_DATA(1:nr:row_all, number_index);   
 Ion_TOF    = ALL_DATA(1:nr:row_all, TOF_index)      *(1e-6); %s  
 Ion_Mass   = ALL_DATA(1:nr:row_all, Mass_index)     *AMU;   
 Ion_Charge = ALL_DATA(1:nr:row_all, Charge_index)   *Qe;   
 Ion_X      = (ALL_DATA(1:nr:row_all, X_index)-Xmid) *(1e-3); %m  
 Ion_Y      = (ALL_DATA(1:nr:row_all, Y_index)-Ymid) *(1e-3); %m  
 Ion_Z      = (ALL_DATA(1:nr:row_all, Z_index)-Zmid) *(1e-3); %m
 Ion_Vx     = ALL_DATA(1:nr:row_all, Vx_index)       *mm_microsec2m_sec;  %m/s 
 Ion_Vy     = ALL_DATA(1:nr:row_all, Vy_index)       *mm_microsec2m_sec;  %m/s
 Ion_Vz     = ALL_DATA(1:nr:row_all, Vz_index)       *mm_microsec2m_sec;  %m/s
 
 Ion_E      = 0.5*Ion_Mass.*((Ion_Vx.^2)+(Ion_Vy.^2)+(Ion_Vz.^2));
 Ion_Vtotal = sqrt((Ion_Vx.^2)+(Ion_Vy.^2)+(Ion_Vz.^2));
 Ion_R      = sqrt((Ion_Y.^2)+(Ion_Z.^2));
 Ion_Vr     = sqrt((Ion_Vy.^2)+(Ion_Vz.^2));
 file_no    = ALL_DATA(1:nr:row_all, file_index);        %
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 n_p = 2;
 n_i = 3;
 n_f = 4; 
 %%% Fragment 1 (ion) initial(i) and final(f) parameters %%%
 clear F1_i_n F1_i_TOF F1_i_Mass F1_i_Charge 
 clear F1_i_X  F1_i_Y F1_i_Z 
 clear F1_i_Vx F1_i_Vy F1_i_Vz
 
 F1_i_n      = ALL_DATA(n_i:nr:row_all, number_index);         
 F1_i_TOF    = (ALL_DATA(n_i:nr:row_all, TOF_index)-Tmid*0) *(1e-6);           
 F1_i_Mass   = ALL_DATA(n_i:nr:row_all, Mass_index)         *AMU;           
 F1_i_Charge = ALL_DATA(n_i:nr:row_all, Charge_index)       *Qe;         
 F1_i_X      = (ALL_DATA(n_i:nr:row_all, X_index)-Xmid)     *(1e-3);            
 F1_i_Y      = (ALL_DATA(n_i:nr:row_all, Y_index)-Ymid)     *(1e-3);             
 F1_i_Z      = (ALL_DATA(n_i:nr:row_all, Z_index)-Zmid)     *(1e-3);         
 F1_i_Vx     = ALL_DATA(n_i:nr:row_all, Vx_index)           *mm_microsec2m_sec;             
 F1_i_Vy     = ALL_DATA(n_i:nr:row_all, Vy_index)           *mm_microsec2m_sec;             
 F1_i_Vz     = ALL_DATA(n_i:nr:row_all, Vz_index)           *mm_microsec2m_sec;             
 
 F1_i_Vx_CM  = F1_i_Vx-Ion_Vx;
 F1_i_Vy_CM  = F1_i_Vy-Ion_Vy;
 F1_i_Vz_CM  = F1_i_Vz-Ion_Vz;
 F1_i_Vo_CM  = sqrt(F1_i_Vx_CM.^2+F1_i_Vy_CM.^2+F1_i_Vz_CM.^2);
 
 F1_i_theta_CM = acos(F1_i_Vx_CM./F1_i_Vo_CM);
 
 
 clear F1_f_n F1_f_TOF F1_f_Mass F1_f_Charge 
 clear F1_f_X  F1_f_Y F1_f_Z 
 clear F1_f_Vx F1_f_Vy F1_f_Vz

 F1_f_n      = ALL_DATA(n_f:nr:row_all, number_index);   
 F1_f_TOF    = (ALL_DATA(n_f:nr:row_all, TOF_index)-Tmid)  *(1e-6);   
 F1_f_Mass   = ALL_DATA(n_f:nr:row_all, Mass_index)        *AMU;   
 F1_f_Charge = ALL_DATA(n_f:nr:row_all, Charge_index)      *Qe;   
 F1_f_X      = (ALL_DATA(n_f:nr:row_all, X_index)-Xmid)    *(1e-3);   
 F1_f_Y      = (ALL_DATA(n_f:nr:row_all, Y_index)-Ymid)    *(1e-3);   
 F1_f_Z      = (ALL_DATA(n_f:nr:row_all, Z_index)-Zmid)    *(1e-3);
 F1_f_Vx     = ALL_DATA(n_f:nr:row_all, Vx_index)          *mm_microsec2m_sec;   
 F1_f_Vy     = ALL_DATA(n_f:nr:row_all, Vy_index)          *mm_microsec2m_sec;   
 F1_f_Vz     = ALL_DATA(n_f:nr:row_all, Vz_index)          *mm_microsec2m_sec;
 
 
 F1_f_R = sqrt(((F1_f_X*(1e-3)).^2)+(F1_f_Z.^2));
 
  
 F1_V_scan    = ALL_DATA (n_p:nr:row_all, 2);
 F1_E_release = ALL_DATA (n_p:nr:row_all, 5);
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 n_p = 6;
 n_i = 7;
 n_f = 8; 
 %%% Fragment 2 initial(i) and final(f) parameters %%%
 clear F2_i_n F2_i_TOF F2_i_Mass F2_i_Charge 
 clear F2_i_X  F2_i_Y F2_i_Z 
 clear F2_i_Vx F2_i_Vy F2_i_Vz
 clear F2_i_theta F2_i_phi_yz F2_i_phi_xz
 F2_i_n      = ALL_DATA(n_i:nr:row_all, number_index);   
 F2_i_TOF    = (ALL_DATA(n_i:nr:row_all, TOF_index)-Tmid)  *(1e-6);   
 F2_i_Mass   = ALL_DATA(n_i:nr:row_all, Mass_index)        *AMU;   
 F2_i_Charge = ALL_DATA(n_i:nr:row_all, Charge_index)      *Qe;   
 F2_i_X      = (ALL_DATA(n_i:nr:row_all, X_index)-Xmid)    *(1e-3);   
 F2_i_Y      = (ALL_DATA(n_i:nr:row_all, Y_index)-Ymid)    *(1e-3);   
 F2_i_Z      = (ALL_DATA(n_i:nr:row_all, Z_index)-Zmid)    *(1e-3);
 F2_i_Vx     = ALL_DATA(n_i:nr:row_all, Vx_index)          *mm_microsec2m_sec;   
 F2_i_Vy     = ALL_DATA(n_i:nr:row_all, Vy_index)          *mm_microsec2m_sec;   
 F2_i_Vz     = ALL_DATA(n_i:nr:row_all, Vz_index)          *mm_microsec2m_sec;
 
 F2_i_theta  = theta_lab(F2_i_Vx,F2_i_Vy,F2_i_Vz); 
 F2_i_phi_yz = phi_lab(F2_i_Vx,F2_i_Vy,F2_i_Vz); % azimuthal angle Y-->Z
 F2_i_phi_xz = phi_lab(F2_i_Vx,F2_i_Vx,F2_i_Vz); % azimuthal angle X-->Z 
  
 F2_i_Vx_CM    = F2_i_Vx-Ion_Vx;
 F2_i_Vy_CM    = F2_i_Vy-Ion_Vy;
 F2_i_Vz_CM    = F2_i_Vz-Ion_Vz;
 F2_i_Vo_CM    = sqrt(F2_i_Vx_CM.^2+F2_i_Vy_CM.^2+F2_i_Vz_CM.^2);
 F2_i_theta_CM = acos(F2_i_Vx_CM./F2_i_Vo_CM);
 
F2_i_E_CM_eV = 0.5*F2_i_Mass.*F2_i_Vo_CM.^2/Qe;
 
 clear F2_f_n F2_f_TOF F2_f_Mass F2_f_Charge 
 clear F2_f_X  F2_f_Y F2_f_Z 
 clear F2_f_Vx F2_f_Vy F2_f_Vz
 F2_f_n      = ALL_DATA(n_f:nr:row_all, number_index);   
 F2_f_TOF    = (ALL_DATA(n_f:nr:row_all, TOF_index)-Tmid)  *(1e-6);   
 F2_f_Mass   = ALL_DATA(n_f:nr:row_all, Mass_index)        *AMU;   
 F2_f_Charge = ALL_DATA(n_f:nr:row_all, Charge_index)      *Qe;   
 F2_f_X      = (ALL_DATA(n_f:nr:row_all, X_index)-Xmid)    *(1e-3);   
 F2_f_Y      = (ALL_DATA(n_f:nr:row_all, Y_index)-Ymid)    *(1e-3);   
 F2_f_Z      = (ALL_DATA(n_f:nr:row_all, Z_index)-Zmid)    *(1e-3);
 F2_f_Vx     = ALL_DATA(n_f:nr:row_all, Vx_index)          *mm_microsec2m_sec;   
 F2_f_Vy     = ALL_DATA(n_f:nr:row_all, Vy_index)          *mm_microsec2m_sec;   
 F2_f_Vz     = ALL_DATA(n_f:nr:row_all, Vz_index)          *mm_microsec2m_sec;
 
 
 
  F2_f_R  = sqrt(((F2_f_X*1e-3).^2)+(F2_f_Z*1e-3.^2));
  F2_f_R2 = (((F2_f_X).^2)+(F2_f_Z.^2)); 
  %F2_f_phi_xz  = phi_lab(F2_f_Y,F2_f_Y,F2_f_Z);
  %F2_f_theta=(F2_f_R*1e3)./(F2_f_TOF*1e9);

  
F2_V_scan    = ALL_DATA (n_p:nr:row_all, 2);
F2_E_release = ALL_DATA (n_p:nr:row_all, 5); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %do_assignment_of_data
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
%% Velocity of ions from start to dissociation

%length = (Ion_X-F2_i_X).*1e-6;
%time_sim = (Ion_TOF-F2_i_TOF)*1e-6;
%v_sim = length./time_sim;

% amu2eV = 9.3149410242e8; % eV/amu
% c = 299792458; %m/s 
% mass =39.9480*amu2eV/c^2; %amuCO
% E_start = 0.1000; %eV
% dE_start  = 0.05; %eV
% v_theory_mean = sqrt(2*E_start/mass);
% v_theory_upper_lim = sqrt(2*(E_start+dE_start)/mass);

%time_theory = length./v_theory_mean;

% figure(1000);clf
%     set(gca,'FontSize',14)
%     hold on; box on
%     plot(v_sim,'.k')
%     plot([0 12500],[v_theory_mean v_theory_mean],'-r','LineWidth',2)
%     plot([0 12500],[v_theory_upper_lim v_theory_upper_lim],'--r','LineWidth',2)
%    % plot([0 12500],[v_theory_lower_lim v_theory_lower_lim],'--r','LineWidth',2)
%     legend('Simulation', 'v(E)','v(E+dE)','interpreter','latex')
%     ylabel('velocity [m/s]','interpreter','latex','FontSize',14)
%     xlabel('ion number','interpreter','latex','FontSize',14)
%     xlim([0 12500])
%     hold off
%     
% figure(1001);clf
%     set(gca,'FontSize',14)
%     hold on; box on
%     %plot(time_sim,'-r','LineWidth',.5)
%     %plot(time_theory,'-k','LineWidth',.5)
%     plot(time_sim-time_theory)
%     legend('Simulation', 'Theory','interpreter','latex')
%     ylabel('Time [s]','interpreter','latex','FontSize',14)
%     xlabel('ion number','interpreter','latex','FontSize',14)
%     xlim([0 12500])
%     hold off




%% Computing gaussian weights

sigma_y = 0.26*1e-3;                             % gaussian width in y %
sigma_z = 0.15*1e-3;                             % gaussian width in z %

z_center = 0;%0.5*1e-3;                         % horizontal center position of peak
y_center = 0;%-1*1e-3;                          % vertical center position of peak
 
z_pos_ion = F1_i_Z;                             % horizontal pos of ions used for weight comp
y_pos_ion = F1_i_Y;                             % vertical pos of ions used for weight comp

w_N = exp(-((z_pos_ion-z_center).^2)./(2*sigma_z^2))...
    .*exp(-((y_pos_ion-y_center).^2)./(2*sigma_y^2));  % total weight of particle %

% Normalizing weights 
 Norm_ion = sum(w_N);                   % norm of weights
 N_ion = length(F1_i_Z);
 w_N = w_N/Norm_ion*N_ion;              % normalized weight so sum(w_N) =N %

%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Figures with the initial distribution %%%%%%%%%%%%%%%%%%
 Y_DET_UP = 1.156%0.139; %0.142
 Y_DET_DW = -0.142%-0.139; %-0.142
 
 cut_ion_detected = (F1_f_Y<Y_DET_DW+0.001).*(F1_f_Y>Y_DET_DW-0.001).*(F1_f_R<=0.02);
 I_ion_detected   = find(cut_ion_detected);
 
 cut_ele_detected = (F2_f_Y<Y_DET_UP+0.001).*(F2_f_Y>Y_DET_UP-0.001).*(F2_f_R<=0.02);
 I_ele_detected   = find(cut_ele_detected);
 
 
 cut_E1 = (F2_i_E_CM_eV>0.1).*(F2_i_E_CM_eV<2.0);
 cut_E2 = (F2_i_E_CM_eV>0.1).*(F2_i_E_CM_eV<2.6);
 cut_E3 = (F2_i_E_CM_eV>0.1).*(F2_i_E_CM_eV<3.6);
 
 cut_ele_E1       = cut_ele_detected.*cut_E1;
 cut_ele_E2       = cut_ele_detected.*cut_E2;
 cut_ele_E3       = cut_ele_detected.*cut_E3;
 I_ele_E1         = find(cut_ele_E1);
 I_ele_E2         = find(cut_ele_E2);
 I_ele_E3         = find(cut_ele_E3);
 
 eff_E1 =sum(cut_ele_E1)/sum(cut_E1)*100;
 eff_E2 =sum(cut_ele_E2)/sum(cut_E2)*100;
 eff_E3 =sum(cut_ele_E3)/sum(cut_E3)*100;
 
  cut_ion_ele_coin = cut_ion_detected.*cut_ele_detected;
  I_ion_ele_coin   = find(cut_ion_ele_coin);
 
%%% some first figures
figure(20);clf
hold on
    set(gca,'FontSize',14)
    title(['eff = ' num2str(eff_E1) ' %'])
    plot(F2_i_X(I_ele_detected)*1e3,F2_i_Z(I_ele_detected)*1e3,'.r')
    plot(F2_f_X(I_ele_detected)*1e3,F2_f_Z(I_ele_detected)*1e3,'.b')
    t=(0:0.1:2*pi);x=40*cos(t); y=40*sin(t);
    plot(x,y,'-')
    xlabel('Longitudinal direction X (mm)','interpreter','latex')
    ylabel('Transverse direction Z (mm)','interpreter','latex')

hold off
box on

tof_bin = (0:1:5000);
hist_F2_TOF = hist((F2_f_TOF(I_ele_detected)-F2_i_TOF(I_ele_detected))*1e9,tof_bin);
hist_F2_TOF_E1 = hist((F2_f_TOF(I_ele_E1)-F2_i_TOF(I_ele_E1))*1e9,tof_bin);
hist_F2_TOF_E2 = hist((F2_f_TOF(I_ele_E2)-F2_i_TOF(I_ele_E2))*1e9,tof_bin);
hist_F2_TOF_E3 = hist((F2_f_TOF(I_ele_E3)-F2_i_TOF(I_ele_E3))*1e9,tof_bin);



A_in_90_91_92 = load('O:\Nat_IFA-xring\CROSSED_BEAMS\ANALYSIS_DATA\O2_april2024\test_save_90_91_92.mat')
%A_in_218 = load('O:\Nat_IFA-xring\CROSSED_BEAMS\ANALYSIS_OF_SIMULATION\MMRUN_218_save.mat')
%A_in_202 = load('O:\Nat_IFA-xring\CROSSED_BEAMS\ANALYSIS_OF_SIMULATION\MMRUN_202_save.mat')
ii_scan =find(A_in_90_91_92.U_scan_1==16.4)
t = (100:1000)*1e-9
    E_ion_eV=3000;
    M_ion_amu = 127;
    Ee_eV = 0.435;
    beta=0;
    V_int=13.8;
    V_drift=15.5;
    L=0.37;
y_model = P_t(E_ion_eV,M_ion_amu,Ee_eV,beta,V_int,V_drift,L,t)
figure(21);clf
hold on
    set(gca,'FontSize',14)
    plot(tof_bin,1*hist_F2_TOF/sum(hist_F2_TOF),'b','linewidth',3)
    %plot(A_in_218.tof_bin,A_in_218.hist_F2_TOF/A_in_218.sum_hist_F2_TOF,'-r','linewidth',3) 
    %plot(A_in_152.tof_bin,A_in_152.hist_F2_TOF/A_in_152.sum_hist_F2_TOF,'-g','linewidth',3)
    %change name for run_ii
    plot(A_in_90_91_92.tof_bin_REF*1e6-75,25*A_in_90_91_92.Hists_U_ESP(ii_scan,:)/sum(A_in_90_91_92.Hists_U_ESP(ii_scan,:)),'g-')
    %plot(tof_bin,hist_F2_TOF_E1,'b')
    %plot(tof_bin,hist_F2_TOF_E2,'g')
    %plot(tof_bin,hist_F2_TOF_E3,'c')
    %plot(t*1e9,y_model/sum(y_model),'-k')
    %plot((F2_f_TOF(I_ele_detected)-F2_i_TOF(I_ele_detected))*1e9,F2_f_Z(I_ele_detected)*1e3,'.b')
    legend('355nm - 0.135eV, Drift=15.5','355nm - Run-90-91-92' ,'Run-0096-97-98 355nm, Drift=16.4', '355nm - Theoretical curve, Drift=21','532nm - 0.867, Drift=20')
    
    xlabel('TOF electron (ns)','interpreter','latex')
    xlim([0 1000])
   
hold off
box on
sum_mikkel = tof_bin,hist_F2_TOF_E1 + tof_bin,hist_F2_TOF_E2 + tof_bin,hist_F2_TOF_E3;
disp(sum_mikkel);



 return
 
 %%
 
 x_bin=(-35:0.4:35)*1e-3;
 y_bin=(-10:0.05:10)*1e-3;
 z_bin=(-10:0.05:10)*1e-3;
 
 % Computing direct and weighted histograms
[hist_direct_all_x, hist_weight_all_x] = uniform_2_gauss(F1_i_X,w_N,x_bin);
[hist_direct_all_y, hist_weight_all_y] = uniform_2_gauss(F1_i_Y,w_N,y_bin);
[hist_direct_all_z, hist_weight_all_z] = uniform_2_gauss(F1_i_Z,w_N,z_bin);

[hist_direct_ion_x, hist_weight_ion_x] = uniform_2_gauss(F1_i_X(I_ion_detected),w_N(I_ion_detected),x_bin);
[hist_direct_ion_y, hist_weight_ion_y] = uniform_2_gauss(F1_i_Y(I_ion_detected),w_N(I_ion_detected),y_bin);
[hist_direct_ion_z, hist_weight_ion_z] = uniform_2_gauss(F1_i_Z(I_ion_detected),w_N(I_ion_detected),z_bin);

[hist_direct_ele_x, hist_weight_ele_x] = uniform_2_gauss(F2_i_X(I_ele_detected),w_N(I_ele_detected),x_bin);
[hist_direct_ele_y, hist_weight_ele_y] = uniform_2_gauss(F2_i_Y(I_ele_detected),w_N(I_ele_detected),y_bin);
[hist_direct_ele_z, hist_weight_ele_z] = uniform_2_gauss(F2_i_Z(I_ele_detected),w_N(I_ele_detected),z_bin);
 
[hist_direct_coin_x, hist_weight_coin_x] = uniform_2_gauss(F2_i_X(I_ion_ele_coin),w_N(I_ion_ele_coin),x_bin);
[hist_direct_coin_y, hist_weight_coin_y] = uniform_2_gauss(F2_i_Y(I_ion_ele_coin),w_N(I_ion_ele_coin),y_bin);
[hist_direct_coin_z, hist_weight_coin_z] = uniform_2_gauss(F2_i_Z(I_ion_ele_coin),w_N(I_ion_ele_coin),z_bin);
 

 
  %%% Numbers for calibration %%%
 N_all_here  = size(F1_i_n,1)
 N_ion_here  = size(I_ion_detected,1)
 N_ele_here  = size(I_ele_detected,1)
 N_coin_here = size(I_ion_ele_coin,1)

 N_all_direct  =sum(hist_direct_all_x)
 N_ion_direct  =sum(hist_direct_ion_x)
 N_ele_direct  =sum(hist_direct_ele_x)
 N_coin_direct =sum(hist_direct_coin_x)

 N_all_weight  =sum(hist_weight_all_x)
 N_ion_weight  =sum(hist_weight_ion_x)
 N_ele_weight  =sum(hist_weight_ele_x)
 N_coin_weight =sum(hist_weight_coin_x)

 
%  L_all_direct  = (max(F1_i_X)-min(F1_i_X))*1e3;      %   
%  L_ion_direct  = (N_ion_direct/N_all_direct)*L_all_direct;  % file 1-20: 31.6 ; file 1-19: %
%  L_ele_direct  = (N_ele_direct/N_all_direct)*L_all_direct;  % file 1_20: 6.9  ; file 1-19:%
%  L_coin_direct = (N_coin_direct/N_all_direct)*L_all_direct; % file 1-20: 6.9  ; file 1-19:%
% 
%  L_all_weight  = (max(F1_i_X)-min(F1_i_X))*1e3;      % 
%  L_ion_weight  = (N_ion_weight/N_all_weight)*L_all_weight;  % file 1-20: 33.7; file 1-19:33.6%
%  L_ele_weight  = (N_ele_weight/N_all_weight)*L_all_weight;  % file 1-20: 5.9 ; file 1-19:5.9%
%  L_coin_weight = (N_coin_weight/N_all_weight)*L_all_weight; % file 1-20: 5.9 ; file 1-19:5.9%
%  
 %%%%%

 
 
 figure(20);clf
 subplot(3,1,1)
 hold on
  %  plot(x_bin*1e3,hist_direct_all_x)
  %  plot(x_bin*1e3,hist_weight_all_x)
    
    plot(x_bin*1e3,hist_direct_ion_x)
    plot(x_bin*1e3,hist_weight_ion_x)
    
    %plot(x_bin*1e3,hist_direct_ele_x,'k')
    plot(x_bin*1e3,hist_weight_ele_x,'k')
    
    plot(x_bin*1e3,hist_direct_coin_x,'g')
    plot(x_bin*1e3,hist_weight_coin_x,'g')
    
 hold off
 box on
 
 subplot(3,1,2)
 hold on
    plot(y_bin*1e3,hist_direct_all_y)
     plot(y_bin*1e3,hist_weight_all_y)
 hold off
 box on
 
 subplot(3,1,3)
 hold on
    plot(z_bin*1e3,hist_direct_all_z)
     plot(z_bin*1e3,hist_weight_all_z)
 hold off
 box on
 
 
 
 
 
 
 
 
 
 
% Distribution of weights
figure(2); clf; 
    hold on; 
    set(gca,'FontSize',14)
    plot3(F1_i_Y*1e3,F1_i_Z*1e3,w_N,'.k')
   % plot3(F1_i_Y(I_ion_detected)*1e3,F1_i_Z(I_ion_detected)*1e3,w_N(I_ion_detected),'.r')
   % plot3(F2_i_Y(I_ele_detected)*1e3,F2_i_Z(I_ele_detected)*1e3,w_N(I_ele_detected),'.b')
    xlabel('Vertical position [mm]','interpreter','latex')
    ylabel('Horizontal position [mm]','interpreter','latex')
    zlabel('Weight','interpreter','latex')
    box on; grid on
    hold off
    
%% figure with starting positions - FIGURE FOR PAPER
c_init = [1 1 1]*0.65%[0 0.4470 0.7410]
c_ion  = [1 0 0]*0.75          %[0.6350 0.0780 0.1840]
c_ele =  [0 0 1]*0.75 %[0.3660 0.6740 0.1880]
c_coin = [1 1 1]*0.5
figure(10);clf
left = 0.2 ; bottom = 0.7 ; width = 0.65; height = 0.25
subplot('position',[left bottom width height])
x_bin_here = (-35:0.33:35)

%ha = tight_subplot(2,1,[.02 .05],[.1 .08],[.13 .02]);
%axes(ha(1));
    hold on
    set(gca,'FontSize',14)
    xticks([])
    xticklabels([-30 -20 -10 0 10 20 30])
    yticks([-2 -1 0 1 2])
    yticklabels([-2 -1 0 1 2])
    h_all = plot(F1_i_X*1e3,F1_i_Y*1e3,'.','color',c_init,'MarkerSize',9);
    h_ion = plot(F1_i_X(I_ion_detected)*1e3,F1_i_Y(I_ion_detected)*1e3,'.','color',c_ion,'MarkerSize',9);
    h_ele = plot(F2_i_X(I_ele_detected)*1e3,F2_i_Y(I_ele_detected)*1e3,'.','color',c_ele,'MarkerSize',9);

    
    
    
    %    h_coin = plot(F1_i_X(I_ion_ele_coin)*1e3,F1_i_Y(I_ion_ele_coin)*1e3,'ok','MarkerSize',2);

    plot([-1 1]*35,[0 0],'--k')
    
    sigma_y=0.26;
    plot([-35 30],[1 1]*sigma_y*2,'--k')
    plot([-35 30],[-1 -1]*sigma_y*2,'--k')
    text(32.5,sigma_y*2, '2$\sigma_{\gamma.y}$','interpreter','latex','Fontsize',14)
    text(32,-sigma_y*2,'-2$\sigma_{\gamma.y}$','interpreter','latex','Fontsize',14)
    
    leg1  = legend([h_all h_ion h_ele],'Initiated','iDET','eDET');
    set(leg1,'interpreter','latex','fontsize',12,'Orientation','horizontal')
    legend boxoff
     
     xlim([-32 40])
     ylim([-2 2.5])
   
    %xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('$y$ (mm)','interpreter','latex')  
%     title('Ionization positions','interpreter','latex','FontSize',18)
    text(-30,2.0,'(a)','interpreter','latex','FontSize',14)
    hold off
    box on

left = 0.2 ; bottom = 0.42; width = 0.65; height = 0.25
subplot('position',[left bottom width height])
    hold on
    set(gca,'FontSize',14)
    xticks([])
    xticklabels([-30 -20 -10 0 10 20 30])
    yticks([-2 -1 0 1 2])
    yticklabels([-2 -1 0 1 2])
    h_all = plot(F1_i_X*1e3,F1_i_Z*1e3,'.','color',c_init,'MarkerSize',9);
    h_ion = plot(F1_i_X(I_ion_detected)*1e3,F1_i_Z(I_ion_detected)*1e3,'.','color',c_ion,'MarkerSize',9);
    h_ele = plot(F2_i_X(I_ele_detected)*1e3,F2_i_Z(I_ele_detected)*1e3,'.','color',c_ele,'MarkerSize',8);
    plot([-1 1]*32,[0 0],'--k')
    
    sigma_z=0.15;
    plot([-35 30],[1 1]*sigma_z*2,'--k')
    plot([-35 30],[-1 -1]*sigma_z*2,'--k')
    text(32.5,sigma_z*2, '2$\sigma_{\gamma.z}$','interpreter','latex','Fontsize',14)
    text(32,-sigma_z*2,'-2$\sigma_{\gamma.z}$','interpreter','latex','Fontsize',14)
    
    
    % h_coin = plot(F1_i_X(I_ion_ele_coin)*1e3,F1_i_Z(I_ion_ele_coin)*1e3,'ok','MarkerSize',2);
    %leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    %set(leg1,'interpreter','latex')
    %xlabel('Longitudinal position $x$ (mm)','interpreter','latex')
    ylabel('$z$ (mm)','interpreter','latex') 
    xlim([-32 40])
    
   
    ylim([-2 2.5])
    
    text(-30,2.0,'(b)','interpreter','latex','FontSize',14)
    box on
    hold off

    left = 0.2 ; bottom = 0.1 ; width = 0.65; height = 0.30
    subplot('position',[left bottom width height])
    hold on
    set(gca,'FontSize',14)
    xticks([-30 -20 -10 0 10 20 30])
    xticklabels([-30 -20 -10 0 10 20 30])
    %yticks([-2 -1 0 1 2])
    %yticklabels([-2 -1 0 1 2])
    %hist_all = hist(F1_i_X*1e3,x_bin_here);
    %S_all = sum(hist_all)
    %plot(x_bin_here,hist_all,'color',[0 0.4470 0.7410],'linewidth',2)
    
%     hist_ion = hist(F1_i_X(I_ion_detected)*1e3,x_bin_here);
%     S_ion = sum(hist_ion)
%     plot(x_bin_here,hist_ion,'color',[0.6350 0.0780 0.1840],'linewidth',2)
%     
%     hist_ele = hist(F1_i_X(I_ele_detected)*1e3,x_bin_here);
%     S_ele = sum(hist_ele)
%     plot(x_bin_here,hist_ele,'color',[0.3660 0.6740 0.1880],'linewidth',2)
    
%     hist_coin = hist(F1_i_X(I_ion_ele_coin)*1e3,x_bin_here,'linewidth',2);
%     S_coin = sum(hist_coin)
%     plot(x_bin_here,hist_coin,'--k')
    
    off_set=0.02
    N_all    = sum(hist_direct_all_x);
    N_weight = sum(hist_weight_all_x);
    
    h1=plot(x_bin*1e3,hist_direct_all_x/N_all+off_set,':','color',c_init,'linewidth',2)
    h2=plot(x_bin*1e3,hist_direct_ion_x/N_all+off_set,'color',c_ion,'linewidth',2)
    h3=plot(x_bin*1e3,hist_direct_ele_x/N_all+off_set,'color',c_ele,'linewidth',2)
    %h4=plot(x_bin*1e3,hist_direct_coin_x/N_all+off_set,'--','color',c_coin,'linewidth',2)
    
    
     
    h1a=plot(x_bin*1e3,hist_weight_all_x/N_weight,':','color',c_init,'linewidth',2)
    h2a=plot(x_bin*1e3,hist_weight_ion_x/N_weight,'-','color',c_ion,'linewidth',2)
    h3a=plot(x_bin*1e3,hist_weight_ele_x/N_weight,'color',c_ele,'linewidth',2)
    %h4a=plot(x_bin*1e3,hist_weight_coin_x/N_weight,'--','color',c_coin,'linewidth',2)
    
    %%% L_ion
        plot([-26.2 -17.4],[1 1]*0.015,'-k')
        plot([-26.2 -25.0],[1 1.05]*0.015,'-k')
        plot([-26.2 -25.0],[1 0.95]*0.015,'-k')
        plot([-18.6 -17.4],[1.05 1]*0.015,'-k')
        plot([-18.6 -17.4],[0.95 1]*0.015,'-k')
        
        plot(-[-26.2 -17.4],[1 1]*0.015,'-k')
        plot(-[-26.2 -25.0],[1 1.05]*0.015,'-k')
        plot(-[-26.2 -25.0],[1 0.95]*0.015,'-k')
        plot(-[-18.6 -17.4],[1.05 1]*0.015,'-k')
        plot(-[-18.6 -17.4],[0.95 1]*0.015,'-k')
        
        plot([-8.6 8.6],[1 1]*0.015,'-k')
        plot([-8.6 -8.6+1.2],[1 1.05]*0.015,'-k')
        plot([-8.6 -8.6+1.2],[1 0.95]*0.015,'-k')
        plot([8.6-1.2 8.6],[1.05 1]*0.015,'-k')
        plot([8.6-1.1 8.6],[0.95 1]*0.015,'-k')
        
    text(30,0.015,'$L_i$','interpreter','latex','fontsize',14) 
    
    %%% L_electron
%         plot([-8.6 8.6],[1 1]*0.016,'-g')
%         plot([-8.6 -8.6+1.2],[1 1.05]*0.016,'-g')
%         plot([-8.6 -8.6+1.2],[1 0.95]*0.016,'-g')
%         plot([8.6-1.2 8.6],[1.05 1]*0.016,'-g')
%         plot([8.6-1.1 8.6],[0.95 1]*0.016,'-g')
%          text(30,0.016,'$L_e$','interpreter','latex','fontsize',14) 
    
    %leg1  = legend([h1 h2 h3 h4],'Initiated','iDET','eDET','iDET $\&$ eDET');
    leg1  = legend([h1 h2 h3],'Initiated','iDET','eDET');
     set(leg1,'interpreter','latex','fontsize',12,'Orientation','horizontal')
     legend boxoff
     
    text(31.5,0.025,['Unif.' newline 'beam'],'interpreter','latex','fontsize',14)
    text(31.5,0.007,['Ph.' newline 'beam'],'interpreter','latex','fontsize',14)
    
   
    
    
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('Normalized intensity','interpreter','latex') 
    xlim([-32 40])
    ylim([0 0.035])
     text(-30,0.031,'(c)','interpreter','latex','FontSize',14)
    hold off
    box on
    
    
%     orient portrait
%           fig_width =6.0;          %5
%           fig_height=8.0;         %6
%            papersize=get(gcf,'Papersize')
%            set(gcf,'PaperType','A4')
%            left=(papersize(1)-fig_width)/2;
%            bottom=(papersize(1)-fig_height)/2+2;
%            set(gcf,'PaperPosition',[left, bottom,fig_width,fig_height])
%     %  print('-depsc2','fig_simulation1.eps'); %color
%     
    
    
    
    mean(F1_i_X(I_ele_detected)*1e3)
    
    %     saveas(gcf,'O:\Nat_IFA-xring\PHOTON_DIAGNOSTIC\SOFIE_MASTERS_THESIS\ionization_pos','epsc')
   
 %%
 clear sig_U sig_Ion
 for jj = 1:length(files)
 jj_file = files(jj)
 cut_ion_detected = (F1_f_Y<Y_DET_DW+0.001).*(F1_f_Y>Y_DET_DW-0.001).*(F1_f_R<=0.02).*(file_no==jj_file);
 I_ion_detected   = find(cut_ion_detected);
 
 cut_ele_detected = (F2_f_Y<Y_DET_UP+0.001).*(F2_f_Y>Y_DET_UP-0.001).*(F2_f_R<=0.02).*(file_no==jj_file);
 I_ele_detected   = find(cut_ele_detected);
 
 cut_ion_ele_coin = cut_ion_detected.*cut_ele_detected;
 I_ion_ele_coin   = find(cut_ion_ele_coin);
 
 
 
 TOF_coin_all = F2_f_TOF(I_ion_ele_coin)- F1_f_TOF(I_ion_ele_coin); %tof_e - tof_i%

 
 tof_bin=(-10:0.03:0)*1e-6;
 [hist_coin_all, hist_coin_weight] = uniform_2_gauss(TOF_coin_all,w_N(I_ion_ele_coin),tof_bin);
 
 
 sig_U(jj) = sum(hist_coin_weight);
 sig_Ion(jj) = length(I_ion_detected);
 end
 
 figure(30);clf
 subplot(3,1,1)
 hold on
    plot(tof_bin*1e6,hist_coin_all)
    plot(tof_bin*1e6,hist_coin_weight)
    plot([-1 -1]*5.68,[0 200],'--k')
    xlim([-8 0])
 hold off
 box on
 
 subplot(3,1,2)
 hold on
    plot(pot,sig_U,'.')
    plot(pot,sig_Ion,'.')
    xlim([-40 0])
 hold off
 box on
 
 subplot(3,1,3)
 hold on
    plot(pot,sig_U./sig_Ion,'.')
    
    xlim([-40 0])
 hold off
 box on
 
 
 
 
    
    
 return 
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(110);clf
     hold on
    set(gca,'FontSize',14)
    h_all = plot(F1_i_X*1e3,F1_i_Z*1e3,'.','color',[0 0.4470 0.7410],'MarkerSize',9);
    h_ion = plot(F1_i_X(I_ion_detected)*1e3,F1_i_Y(I_ion_detected)*1e3,'.','color',[0.6350 0.0780 0.1840],'MarkerSize',9);
    h_ele = plot(F2_i_X(I_ele_detected)*1e3,F2_i_Y(I_ele_detected)*1e3,'.','color',[0.3660 0.6740 0.1880],'MarkerSize',9);
   % h_coin = plot(F1_i_X(I_ion_ele_coin)*1e3,F1_i_Y(I_ion_ele_coin)*1e3,'ok','MarkerSize',2);
%     leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
%     set(leg1,'interpreter','latex','fontsize',12)
    xlim([-35 35])
    xticks([-30 -20 -10 0 10 20 30])
    xticklabels([-30 -20 -10 0 10 20 30])
    ylim([-2.5 2.5])
    yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
    yticklabels([-1.5 -1 -0.5 0 0.5 1 1.5])
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Vertical position [mm]','interpreter','latex')  
%     title('Ionization positions','interpreter','latex','FontSize',18)
    leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons','Position',[400 215 0.25 0.1]);
    set(leg1,'interpreter','latex','fontsize',14)

    hold off
    box on
 return   
   %%
 figure(11); clf
 hold on
    set(gca,'FontSize',14)
    h_all = stairs(x_bin*1e3,hist_x_ion_all,'-b');
    h_ion = stairs(x_bin*1e3,hist_x_ion_detected,'-r');
    h_trap = stairs(x_bin*1e3,hist_x_ele_detected,'--g','Linewidth',1.5);
    leg1  = legend([h_all h_ion h_trap],'All','Detected ions','Detected electrons');
    set(leg1,'interpreter','latex')
    %legend boxoff
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Intensity [counts]','interpreter','latex')
     
     xlim([-40 40])
     %ylim([-5 50])
 %hold off
 box on
 print(gcf,'vertical_distribution.png')

 figure(12); clf
 subplot(2,1,1)
 hold on
    set(gca,'FontSize',14)
    h_all = plot(F1_i_Y*1e3,F1_i_Z*1e3,'.b');
    h_ion = plot(F1_i_Z(I_ion_detected)*1e3,F1_i_Y(I_ion_detected)*1e3,'.r');
    
    h_ele = plot(F2_i_Z(I_ele_detected)*1e3,F2_i_Y(I_ele_detected)*1e3,'og','MarkerSize',2);
    
    leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    set(leg1,'interpreter','latex')
    legend boxoff
    xlim([-10 10])
    ylim([-5 5])
    xlabel('Z position [mm]','interpreter','latex')
    ylabel('Y position [mm]','interpreter','latex')   
%hold off
 box on
 
 subplot(2,1,2)
 hold on
    set(gca,'FontSize',14)
    h_all = stairs(y_bin*1e3,hist_y_ion_all,'-b');
    h_ion = stairs(y_bin*1e3,hist_y_ion_detected,'-r');
    h_ele = stairs(y_bin*1e3,hist_y_ele_detected,'--g','Linewidth',1.5);
    leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    set(leg1,'interpreter','latex')
    legend boxoff
     xlabel('Longitudinal position [mm]','interpreter','latex')
     xlabel('Intensity [counts]','interpreter','latex')
     
     xlim([-5 5])
     %ylim([-5 50])
 %hold off
 box on

%% Final positions

 figure(1020);clf
 subplot(2,1,1)
 hold on
    set(gca,'FontSize',14)
    h_all = plot(F2_f_X*1e3,F2_f_Y*1e3,'.');
    %h_ion = plot(F1_f_X(I_ion_detected)*1e3,F1_f_Y(I_ion_detected)*1e3,'.r');
    h_ele = plot(F2_f_X(I_ele_detected)*1e3,F2_f_Y(I_ele_detected)*1e3,'og','MarkerSize',2);
   % h_coin = plot(F1_i_X(I_ion_ele_coin)*1e3,F1_i_Y(I_ion_ele_coin)*1e3,'ok','MarkerSize',2);
    %leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    %set(leg1,'interpreter','latex')
    %xlim([-35 35])
    %ylim([-3 3])
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Vertical position [mm]','interpreter','latex')   
    hold off
    box on
 subplot(2,1,2)
    hold on
    set(gca,'FontSize',14)
    h_all = plot(F2_f_X*1e3,F2_f_Z*1e3,'.');
    %h_ion = plot(F1_f_X(I_ion_detected)*1e3,F1_f_Z(I_ion_detected)*1e3,'.r');
    h_ele = plot(F2_f_X(I_ele_detected)*1e3,F2_f_Z(I_ele_detected)*1e3,'og','MarkerSize',2);
   % h_coin = plot(F1_i_X(I_ion_ele_coin)*1e3,F1_i_Z(I_ion_ele_coin)*1e3,'ok','MarkerSize',2);
    %leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    %set(leg1,'interpreter','latex')
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Horizontal position [mm]','interpreter','latex') 
    %xlim([-35 35])
    %ylim([-3 3])
    box on
    hold off
%% Check for trapped particles

% Electrons
mean_TOF_ele_det = mean(F2_f_TOF(I_ele_detected)-F2_i_TOF((I_ele_detected)));
trapped_ele = find((F2_f_TOF-F2_i_TOF)>=2*mean_TOF_ele_det); %(F2_f_Y<0.0015)
num_of_ele_trapped = size(trapped_ele,1)

% Ions
mean_TOF_ion_det = mean(F1_f_TOF(I_ion_detected)-F1_i_TOF((I_ion_detected)));
trapped_ion = find((F1_f_TOF-F1_i_TOF)>=2*mean_TOF_ion_det); 
num_of_ions_trapped = size(trapped_ion,1)

% Plots
figure(1004);clf
 subplot(2,1,1)
 hold on
    set(gca,'FontSize',14)
    h_all = plot(F2_i_X*1e3,F2_i_Z*1e3,'.');
    h_ele = plot(F2_i_X(I_ele_detected)*1e3,F2_i_Y(I_ele_detected)*1e3,'ok','MarkerSize',2);
    h_coin = plot(F2_i_X(I_ion_ele_coin)*1e3,F2_i_Y(I_ion_ele_coin)*1e3,'or','MarkerSize',2);
    h_trap = plot(F2_i_X(trapped_ele)*1e3,F2_i_Y(trapped_ele)*1e3,'og','MarkerSize',2);
    leg1  = legend([h_all h_ele h_coin h_trap],'All','Detected electrons','Coincidences','Trapped electrons');
    set(leg1,'interpreter','latex')
    xlim([-35 35])
    ylim([-3 3])
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Vertical position [mm]','interpreter','latex')
    title('Initial positions of electrons','interpreter','latex','FontSize',16)
    hold off
    box on
 subplot(2,1,2)
    hold on
    set(gca,'FontSize',14)
    h_all = plot(F2_i_X*1e3,F2_i_Z*1e3,'.');
    h_ele = plot(F2_i_X(I_ele_detected)*1e3,F2_i_Z(I_ele_detected)*1e3,'ok','MarkerSize',2);
    h_coin = plot(F2_i_X(I_ion_ele_coin)*1e3,F2_i_Z(I_ion_ele_coin)*1e3,'or','MarkerSize',2);
    h_trap = plot(F2_i_X(trapped_ele)*1e3,F2_i_Z(trapped_ele)*1e3,'og','MarkerSize',2);
    leg1  = legend([h_all h_ele h_coin h_trap],'All','Detected electrons','Coincidences','Trapped electrons');
    set(leg1,'interpreter','latex')
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Horizontal position [mm]','interpreter','latex') 
    xlim([-35 35])
    ylim([-3 3])
    box on
    hold off

% figure(1002);clf
%     set(gca,'FontSize',14)
%     hold on; box on
%     plot(F2_f_TOF-F2_i_TOF,'-r','LineWidth',.5)
%     plot(I_ele_detected,F2_f_TOF(I_ele_detected)-F2_i_TOF(I_ele_detected),'-k','LineWidth',.5)
%     plot(trapped_ele,F2_f_TOF(trapped_ele)-F2_i_TOF(trapped_ele),'ob','LineWidth',.5)
%     legend('All', 'Detected','Trapped')
%     set(legend,'interpreter','latex')
%     ylabel('TOF [s]','interpreter','latex','FontSize',14)
%     xlabel('Electron number','interpreter','latex','FontSize',14)
%     xlim([0 12500])
%     hold off
%     
% figure(1003);clf
%     set(gca,'FontSize',14)
%     hold on; box on
%     plot(F2_f_Y*1e3,'-r','LineWidth',.5)
%     plot(I_ele_detected,F2_f_Y(I_ele_detected)*1e3,'-k','LineWidth',.5)
%     plot(trapped_ele,F2_f_Y(trapped_ele)*1e3,'ob','LineWidth',.5)
%     legend('All', 'Detected','Trapped')
%     set(legend,'interpreter','latex')
%     plot([0 12500],[1.5 1.5],'-g')
%     plot([0 12500],[-1.5 -1.5],'-g')
%     ylabel('Final position (y) [mm]','interpreter','latex','FontSize',14)
%     xlabel('Electron number','interpreter','latex','FontSize',14)
%     xlim([0 12500])
%     hold off   

%% conversion from square to gaussian beam - initial positions

% Computing direct and weighted histograms
[hist_direct_all_x, hist_weight_all_x] = uniform_2_gauss(F1_i_X,w_N,x_bin);
[hist_direct_all_y, hist_weight_all_y] = uniform_2_gauss(F1_i_Y,w_N,x_bin);
[hist_direct_all_z, hist_weight_all_z] = uniform_2_gauss(F1_i_Z,w_N,x_bin);

[hist_direct_ion_z, hist_weight_ion_z] = uniform_2_gauss(F1_i_Z(I_ion_detected),w_N(I_ion_detected),x_bin);
[hist_direct_ion_y, hist_weight_ion_y] = uniform_2_gauss(F1_i_Y(I_ion_detected),w_N(I_ion_detected),x_bin);
[hist_direct_ion_x, hist_weight_ion_x] = uniform_2_gauss(F1_i_X(I_ion_detected),w_N(I_ion_detected),x_bin);

[hist_direct_ele_z, hist_weight_ele_z] = uniform_2_gauss(F2_i_Z(I_ele_detected),w_N(I_ele_detected),x_bin);
[hist_direct_ele_y, hist_weight_ele_y] = uniform_2_gauss(F2_i_Y(I_ele_detected),w_N(I_ele_detected),x_bin);
[hist_direct_ele_x, hist_weight_ele_x] = uniform_2_gauss(F2_i_X(I_ele_detected),w_N(I_ele_detected),x_bin);
%%
int = linspace(-2e-3,2e-3,100);
gauss = @(x,center,sigma) exp(-((x-center).^2)./(2*sigma^2));

figure(30);clf
ha = tight_subplot(2,2,[.02 .08],[.1 .01],[.13 .02]);
axes(ha(1));
    hold on
    set(gca,'FontSize',14)
    stairs(x_bin*1e3,hist_direct_ion_z,'-','color',[0.6350 0.0780 0.1840],'linewidth',1.5)
    stairs(x_bin*1e3,hist_weight_ion_z,'-','color',[0 0.4470 0.7410],'linewidth',1.5)
    plot(int*1e3,175*gauss(int,z_center,sigma_z),'--k','linewidth',1.5)
%     xlabel('z position [mm]','interpreter','latex','FontSize',14)
    ylabel('Intensity [counts]','interpreter','latex','FontSize',14,'position',[-2.75 1])
    box on
    %title('Ion','interpreter','latex','FontSize',16)
    txt = strcat('Ion');   
    text(-1.7,220,txt,'interpreter','latex','FontSize',16)
    xlim([-2 2])
    xticks([-1 0 1])
    ylim([0 250])
    yticks([50 100 150 200])
    yticklabels([50 100 150 200])
    hold off
axes(ha(2));
    hold on
    set(gca,'FontSize',14)
    stairs(x_bin*1e3,hist_direct_ion_y,'-','color',[0.6350 0.0780 0.1840],'linewidth',1.5)
    stairs(x_bin*1e3,hist_weight_ion_y,'-','color',[0 0.4470 0.7410],'linewidth',1.5)
    plot(int*1e3,175*gauss(int,y_center,sigma_y),'--k','linewidth',1.5)
%     xlabel('y position [mm]','interpreter','latex','FontSize',14)
%     ylabel('$N_{ion}$ [counts]','interpreter','latex','FontSize',14)
    box on
%     title('Ion','interpreter','latex','FontSize',16)
    txt = strcat('Ion');   
    text(-1.7,220,txt,'interpreter','latex','FontSize',16)
    xlim([-2 2])
    xticks([-1 0 1])
    ylim([0 250])
    yticks([50 100 150 200])
    yticklabels([50 100 150 200])
    hold off
axes(ha(3));
    hold on
    set(gca,'FontSize',14)
    stairs(x_bin*1e3,hist_direct_ele_z,'-','color',[0.6350 0.0780 0.1840],'linewidth',1.5)
    stairs(x_bin*1e3,hist_weight_ele_z,'-','color',[0 0.4470 0.7410],'linewidth',1.5)
    plot(int*1e3,40*gauss(int,z_center,sigma_z),'--k','linewidth',1.5)
    xlabel('Horizontal position [mm]','interpreter','latex','FontSize',14)
%     ylabel('$N_{ele}$ [counts]','interpreter','latex','FontSize',14)
    box on
%     title('Electron','interpreter','latex','FontSize',16)
    txt = strcat('Electron');   
    text(-1.7,60,txt,'interpreter','latex','FontSize',16)
    xlim([-2 2])
    xticks([-1 0 1])
    xticklabels([-1 0 1])
    ylim([0 67])
    yticks([20 40 60])
    yticklabels([20 40 60])
    hold off
axes(ha(4));
    hold on
    set(gca,'FontSize',14)
    stairs(x_bin*1e3,hist_direct_ele_y,'-','color',[0.6350 0.0780 0.1840],'linewidth',1.5)
    stairs(x_bin*1e3,hist_weight_ele_y,'-','color',[0 0.4470 0.7410],'linewidth',1.5)
    plot(int*1e3,75*gauss(int,z_center,sigma_z),'--k','linewidth',1.5)
    xlabel('Vertical position [mm]','interpreter','latex','FontSize',14)
%     ylabel('$N_{ele}$ [counts]','interpreter','latex','FontSize',14)
    box on
    %title('Electron','interpreter','latex','FontSize',16)
    txt = strcat('Electron');   
    text(-1.7,86,txt,'interpreter','latex','FontSize',16)
    xlim([-2 2])
    xticks([-1 0 1])
    xticklabels([-1 0 1])
    ylim([0 95])
    yticks([20 40 60 80])
    yticklabels([20 40 60 80])
    legend('Direct','Weighted','Gaussian')
    set(legend,'interpreter','latex','fontsize',14)
    legend boxoff    
    hold off
    
%     saveas(gcf,'O:\Nat_IFA-xring\PHOTON_DIAGNOSTIC\SOFIE_MASTERS_THESIS\uniform_2_gaussian_beam','epsc')
%%
figure(31); clf;
subplot(2,1,1)
    hold on
    set(gca,'FontSize',14)
    h_all = stairs(x_bin*1e3,hist_direct_all_x,'-b');
    h_ion = stairs(x_bin*1e3,hist_direct_ion_x,'-r');
    h_ele = stairs(x_bin*1e3,hist_direct_ele_x,'--g','Linewidth',1.5);
    leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    set(leg1,'interpreter','latex')
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Intensity [counts]','interpreter','latex')  
    title('Direct histogram')
    xlim([-35 35])
    %ylim([-5 50])
    hold off
    box on
subplot(2,1,2)
    hold on
    set(gca,'FontSize',14)
    h_all = stairs(x_bin*1e3,hist_weight_all_x,'-b');
    h_ion = stairs(x_bin*1e3,hist_weight_ion_x,'-r');
    h_ele = stairs(x_bin*1e3,hist_weight_ele_x,'--g','Linewidth',1.5);
    leg1  = legend([h_all h_ion h_ele],'All','Detected ions','Detected electrons');
    set(leg1,'interpreter','latex')
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Intensity [counts]','interpreter','latex')  
    title('Weighted histogram')
    xlim([-35 35])
    %ylim([-5 50])
    hold off
    box on
% 
%     %% 3D fig of staring positions
% figure(32); clf;
%     hold on; grid on
%     set(gca,'FontSize',14)
%     plot3(F1_i_X.*1e3,F1_i_Y.*1e3,F1_i_Z.*1e3,'.b')
%     plot3(F1_i_X(I_ion_detected).*1e3,F1_i_Y(I_ion_detected).*1e3,F1_i_Z(I_ion_detected).*1e3,'or')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     hold off
% 

%% 3D hist of starting positions (not weighted)
 x_bin=(-100:0.05:100)*1e-3;

hist_ion_i_YZ = [F1_i_Y F1_i_Z];
hist_ele_i_YZ = [F2_i_Y F2_i_Z];

% figure(33); clf;
% subplot(2,1,1)
%     hold on;
%     set(gca,'FontSize',14)
%     hist3(hist_ion_i_YZ,'Ctrs',{-0.002:0.0001:0.002 -0.002:0.0001:0.002},'CdataMode','auto')
%     xlabel('Y')
%     ylabel('Z')
%     xlim([-2.1 2.1]*1e-3)
%     ylim([-2.1 2.1]*1e-3)
%     colorbar
%     hold off
% subplot(2,1,2)
%     hold on;
%     set(gca,'FontSize',14)
%     hist3(hist_ele_i_YZ,'Ctrs',{-0.002:0.0001:0.002 -0.002:0.0001:0.002},'CdataMode','auto')
%     xlabel('Y')
%     ylabel('Z')
%     xlim([-2.1 2.1]*1e-3)
%     ylim([-2.1 2.1]*1e-3)
%     colorbar
%     hold off   

%% Effective lengths

% ---------------------------------- Electrons ----------------------------------
min_ele_det = min(F2_i_X(I_ele_detected));                   % minimum stating pos (X) for detected electrons
max_ele_det = max(F2_i_X(I_ele_detected));                   % maximum stating pos (X) for detected electrons

L_ele_det = abs(min_ele_det)+max_ele_det;                    % length of starting region
N_ele_det = size(I_ele_detected,1);                          % number of detected electrons

cut_ele_in_L = (F2_i_X<max_ele_det).*(F2_i_X>min_ele_det);   % cut all electrons within starting region
I_ele_in_L   = find(cut_ele_in_L);                           % indicies for electrons starting region
N_ele_in_L   = size(I_ele_in_L,1);                           % number of electrons starting in region

[hist_direct_ele_in_L_x, hist_weight_ele_in_L_x] = uniform_2_gauss(F2_i_X(I_ele_in_L),w_N(I_ele_in_L),x_bin);   % hist all electrons in L

% Effective length = #(detected ele)/#(all ele in L)*L
L_eff_ele = hist_weight_ele_x/hist_weight_ele_in_L_x*L_ele_det %N_ele_det/N_ele_in_L*L_ele_det   % effective length (m)


% ------------------------------------ Ions ------------------------------------
min_ion_det = min(F1_i_X(I_ion_detected));                      % minimum stating pos (X) for detected ions
max_ion_det = max(F1_i_X(I_ion_detected));                      % maximum stating pos (X) for detected ions

L_ion_det = abs(min_ion_det)+max_ion_det;                       % length of starting region
N_ion_det = size(I_ion_detected,1);                             % number of detected ions

cut_ion_in_L = (F1_i_X<max_ion_det).*(F1_i_X>min_ion_det);      % cut all ions within starting region
I_ion_in_L   = find(cut_ion_in_L);                              % indicies for ions starting region
N_ion_in_L   = size(I_ion_in_L,1);                              % number of ions starting in region

[hist_direct_ion_in_L_x, hist_weight_ion_in_L_x] = uniform_2_gauss(F1_i_X(I_ion_in_L),w_N(I_ion_in_L),x_bin);  % hist all ions in L

L_eff_ion = hist_weight_ion_x/hist_weight_ion_in_L_x*L_ion_det  % effective length (m)    
    

% -------------------------------- Coincidences -------------------------------- 
% Using ion starting positions but these are the same as for the electrons

min_coin_det = min(F1_i_X(I_ion_ele_coin));                     % minimum stating pos (X) for detected coincidence ions
max_coin_det = max(F1_i_X(I_ion_ele_coin));                     % maximum stating pos (X) for detected coincidence ions

L_coin_det = abs(min_coin_det)+max_coin_det;                    % length of starting region
N_coin_det = size(I_ion_ele_coin,1);                            % number of detected coincidences

cut_coin_in_L = (F1_i_X<max_coin_det).*(F1_i_X>min_coin_det);   % cut all ions within starting region
I_coin_in_L   = find(cut_coin_in_L);                            % indicies for ions starting region
N_coin_in_L   = size(I_coin_in_L,1);                            % number of ions starting in region

[hist_direct_coin_in_L_x, hist_weight_coin_in_L_x] = uniform_2_gauss(F1_i_X(I_coin_in_L),w_N(I_coin_in_L),x_bin);  % hist coincidences in L
[hist_direct_coin_x, hist_weight_coin_x] = uniform_2_gauss(F1_i_X(I_ion_ele_coin),w_N(I_ion_ele_coin),x_bin);      % hist detected coincidences

L_eff_coin = hist_weight_coin_x/hist_weight_coin_in_L_x*L_coin_det  % effective length (m)


 %% TOF for ions and electrons separately
 TOF_elec_pure = F2_f_TOF(I_ele_detected)-F2_i_TOF(I_ele_detected);
 TOF_ion_pure = F1_f_TOF(I_ion_detected)-F1_i_TOF(I_ion_detected);

 TOF_elec_coin = F2_f_TOF(I_ion_ele_coin)-F2_i_TOF(I_ion_ele_coin);
 TOF_ion_coin  = F1_f_TOF(I_ion_ele_coin)-F1_i_TOF(I_ion_ele_coin);
 
 TOF_bin_pure =(-10:0.01:10)*1e-6;
 
 h_ion  = hist(TOF_ion_pure,TOF_bin_pure);
 h_elec = hist(TOF_elec_pure*100,TOF_bin_pure);

 h_ion_coin  = hist(TOF_ion_coin,TOF_bin_pure);
 h_elec_coin = hist(TOF_elec_coin*100,TOF_bin_pure);
 
figure(100);clf
 hold on
     h1=stairs(TOF_bin_pure*1e6,h_ion,'g');
     h2=stairs(TOF_bin_pure*1e6,h_elec,'b');
     h3=stairs(TOF_bin_pure*1e6,h_ion_coin,'k');
     h4=stairs(TOF_bin_pure*1e6,h_elec_coin,'r');
     legend('Ion pure','Electron pure (scaled)','Ion coin','Electron coin (scaled)')
     xlim([0 10])
 hold off
 
 %% detector images
 t   = (0:0.1:2.1*pi);
 x_d = 20*cos(t);
 y_d = 20*sin(t);
 figure(11);clf
 subplot(2,1,1)
 title('Ion detector')
 hold on
  set(gca,'FontSize',14)
    %h_all = plot(F1_f_X*1e3,F1_f_Z*1e3,'.')
    h_ion  = plot(F1_f_X(I_ion_detected)*1e3,F1_f_Z(I_ion_detected)*1e3,'.b');
    h_coin = plot(F1_f_X(I_ion_ele_coin)*1e3,F1_f_Z(I_ion_ele_coin)*1e3,'or');
    %h_ele = plot(F2_f_Y(I_ele_detected)*1e3,F2_f_Z(I_ele_detected)*1e3,'og','MarkerSize',10)
    
    plot(x_d,y_d,'-k','Linewidth',2)
    leg1 = legend([h_ion h_coin],'Detected ions','Coincidences');
    set(leg1,'interpreter','latex')
    legend boxoff
    
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Transverse position [mm]','interpreter','latex')
    
    xlim([-25 25])
    ylim([-25 25])
 %hold off
 box on
 
 
 subplot(2,1,2)
 hold on
 title('Electron detector')
  set(gca,'FontSize',14)
    %h_all = plot(F2_f_X*1e3,F2_f_Z*1e3,'.')
    
    h_ele  = plot(F2_f_X(I_ele_detected)*1e3,F2_f_Z(I_ele_detected)*1e3,'.b');
    h_coin = plot(F2_f_X(I_ion_ele_coin)*1e3,F2_f_Z(I_ion_ele_coin)*1e3,'or');
    
    plot(x_d,y_d,'-k','Linewidth',2)
    
    leg1   = legend([h_ele h_coin],'Detected elec.','Coincidences');
    set(leg1,'interpreter','latex')
    legend boxoff
    xlim([-25 25])
    ylim([-25 25])
    xlabel('Longitudinal position [mm]','interpreter','latex')
    ylabel('Transverse position [mm]','interpreter','latex')
 %hold off
 box on
 
 %% Time of flight
 TOF_bin  = (-10:0.025:10)*1e-6; %(-10:0.025:10)*1e-6
 
 TOF_ion_ele_coin = F2_f_TOF(I_ion_ele_coin)-F1_f_TOF(I_ion_ele_coin);
 
 TOF_ion_coin = F1_f_TOF(I_ion_ele_coin);
 TOF_ele_coin = F2_f_TOF(I_ion_ele_coin);
 
 mean_TOF_ion = mean(TOF_ion_coin)-mean(F1_i_TOF(I_ion_ele_coin));
 mean_TOF_ele = mean(TOF_ele_coin)-mean(F2_i_TOF(I_ion_ele_coin));
 
 hist_TOF_ion_elec_coin=hist(TOF_ion_ele_coin,TOF_bin);
  
 mean_TOF_coinc = mean(TOF_ion_ele_coin)
 TOF_std  = std(TOF_ion_ele_coin)

 % Compute TOF spectrum with gaussian weights
 [hist_direct_ele_TOF, hist_weight_ele_TOF]   = uniform_2_gauss(F2_f_TOF(I_ele_detected)-F2_i_TOF(I_ele_detected),w_N(I_ele_detected),TOF_bin);
 [hist_direct_ion_TOF, hist_weight_ion_TOF]   = uniform_2_gauss(F1_f_TOF(I_ion_detected)-F1_i_TOF(I_ion_detected),w_N(I_ion_detected),TOF_bin);
 [hist_direct_coin_TOF, hist_weight_coin_TOF] = uniform_2_gauss(F2_f_TOF(I_ion_ele_coin)-F1_f_TOF(I_ion_ele_coin),w_N(I_ion_ele_coin),TOF_bin);

 figure(20);clf %Direct and weighted TOF Ar peak
 hold on
 set(gca,'FontSize',14)
    stairs(TOF_bin*1e6,hist_TOF_ion_elec_coin,'b')
    stairs(TOF_bin*1e6,hist_weight_coin_TOF,'r')
    plot([1 1]*mean_TOF_coinc*1e6,[0 600],'--b')
    legend('Direct','Weighted','Direct mean')
    xlabel('TOF$_{electron}$ - TOF$_{ion}$','interpreter','latex')
    ylabel('Intensity [counts]','interpreter','latex')
    %xlim([-10 -2])
 hold off
 box on
 
 sum_hist = sum(hist_weight_coin_TOF)
return 
 figure(21);clf %Direct and weighted electron/ion signals
 hold on
    stairs(TOF_bin*1e6,hist_direct_ele_TOF)
    stairs(TOF_bin*1e6,hist_weight_ele_TOF)
    stairs(TOF_bin*1e6,hist_direct_ion_TOF)
    stairs(TOF_bin*1e6,hist_weight_ion_TOF)
    legend('Electron direct','Electron weighted','Ion pure','Ion weighted')
    %xlim([-1 10])
 hold off

%filename = strcat('TOF_file_232_to_236_1_up.mat');
%save(filename,'hist_weight_coin_TOF','hist_weight_ion_TOF','hist_weight_ele_TOF')

%% Compare sim vs. data
compare_sim_vs_data = 1;
if compare_sim_vs_data==1
    filename = "coin_spectrum_data_Run_059.txt";
    filepath = "O:\Nat_IFA-xring\PHOTON_DIAGNOSTIC\ANALYSIS_of_DATA/Spectrometer/";
    fullpath = filepath+filename;

    data         = load(fullpath);
    tof_coin_bin = data(:,1); 
    signal       = data(:,2);
    background   = data(:,3);
    
    old_sim1      = load('TOF_file_144_to_148_Y0_0_Z0_0.txt');
    old_sim2      = load('TOF_file_182_to_186.txt');
    old_sim3      = load('TOF_file_192_to_196.txt');
    
%      x=(-10:0.055:10).*1e6;
%      y=hist_weight_coin_TOF;
%      fit_sim = fit(x.',y.','gauss2')

    figure(22); clf;
        hold on; box on;
        set(gca,'FontSize',14)
        plot(tof_coin_bin*1e3,signal,'r')
        plot(tof_coin_bin.*1e3,background,'b')
        stairs(TOF_bin*1e12,old_sim1/150+mean(background),'g')
        stairs(TOF_bin*1e12,old_sim2/200+mean(background),'m')
        stairs(TOF_bin*1e12,old_sim3/70+mean(background),'y')
        stairs(TOF_bin*1e12,hist_weight_coin_TOF/100+mean(background),'k')
        legend('Signal 59','Background','Simulation 144-148','Simulation 182-186','Simulation 192-196','Simulation 197-201','Location','northwest')
        xlabel('T$_e$ - T$_{ion}$','interpreter','latex','FontSize',14)
        ylabel('Counts/event','interpreter','latex','FontSize',14)
        xlim([min(tof_coin_bin*1000) max(tof_coin_bin*1000)])
        hold off
end

%%

 
 return
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 clear Xbin Ybin Zbin YZ_hist_ion
 Xbin=(-100:100)*0.01;
 Ybin=(-100:200)*0.01-1;
 Zbin=(-100:200)*0.01-1;
 YZ_hist_ion=histogram_2D_new(F1_i_Z'*(1e3),F1_i_Y'*(1e3),Zbin,Ybin);   
 figure(1);clf
 subplot(2,1,1)
    colormap('default');
    hold on
    pcolor(Zbin,Ybin,YZ_hist_ion);
    colorbar
    shading flat
    xlabel('Ions horizontal coordinate Z(mm)','FontSize',10)
    ylabel('Ions Vertical coordinate Y(mm)','FontSize',10)
    hold off
    box on
 return
 subplot(2,2,2)
    hold on
    hist(Ion_Y*(1e3),Ybin)
    xlabel('Ions Vertical coordinate Y(mm)','FontSize',10)
    axis([-1 1 0 100])
    hold off
    box on
    
 subplot(2,2,3)
    hold on
    hist(Ion_Z*(1e3),Zbin)
    xlabel('Ions horizontal coordinate Z(mm)','FontSize',10)
    axis([-1 1 0 100])
    hold off
    box on
    
  subplot(2,2,4)
    hold on
    %plot(F2_i_X*(1e3),'.')
    hist(Ion_X*(1e3) ,Xbin)
    %hist(F1_i_X,1000)
    xlabel('Ions longitudinal coordinate X(mm)','FontSize',10)
    %axis([-1 1 0 100])
    hold off
    box on   
    
   figure(2);clf
   hold on
   hist(F2_i_theta_CM)
   hold off
   
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% plot of TOF for DET 1
%do_HEX=0
if (do_HEX==1)
det1_Ypos_min=-0.140;
det1_Ypos_max= 0.140;

det1_Zpos_min=-0.140;
det1_Zpos_max= 0.140;

det1_Xpos_min=0.2625;
det1_Xpos_max=0.2635;


 
cut_det1_F1=   (F1_f_Y>=det1_Ypos_min).*(F1_f_Y<=det1_Ypos_max)...
             .*(F1_f_Z>=det1_Zpos_min).*(F1_f_Z<=det1_Zpos_max)...
             .*(F1_f_X>=det1_Xpos_min).*(F1_f_X<=det1_Xpos_max);

alfa=tan(60/180*pi);

cut_star_F2= (F2_f_Y*1e3>-10).*(F2_f_Y*1e3<20+alfa*F2_f_Z*1e3).*(F2_f_Y*1e3<20-alfa*F2_f_Z*1e3)...
            +(F2_f_Y*1e3<10).*(F2_f_Y*1e3>-20+alfa*F2_f_Z*1e3).*(F2_f_Y*1e3>-20-alfa*F2_f_Z*1e3)

 cut_det1_F2=  (F2_f_Y>=det1_Ypos_min).*(F2_f_Y<=det1_Ypos_max)...
            .*(F2_f_Z>=det1_Zpos_min).*(F2_f_Z<=det1_Zpos_max)...
            .*(F2_f_X>=det1_Xpos_min).*(F2_f_X<=det1_Xpos_max)...
            .*(F2_E_release<40)...
            .*(cut_star_F2==0);
 

      
Loss_F1=sum(cut_det1_F1==0)/(sum(cut_det1_F1==1)+sum(cut_det1_F1==0));
Loss_F2=sum(cut_det1_F2==0)/(sum(cut_det1_F2==1)+sum(cut_det1_F2==0));


I_det1_F1=find(cut_det1_F1);
I_det1_F2=find(cut_det1_F2);



 mean_TOF1=mean(F2_i_TOF)*1e9
 mean_TOF2=mean(F2_i_TOF)*1e9

 
figure(5);clf
    subplot(4,1,1)
        hold on
        x=(0:5:12000)
        a1=hist(F1_f_TOF(I_det1_F1)*1e9-mean_TOF2,x);
        a2=hist(F2_f_TOF(I_det1_F2)*1e9-mean_TOF2,x);
  
        plot(x,a1,'r','LineWidth',2)
        plot(x,a2,'b','LineWidth',2)
        %plot(x,a2,'b','LineWidth',3)
        xlim([0000 3000])
        xlabel('TOF (ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        %legend('0.87 eV')
        %legend('1 eV','2.3 eV')
    hold off
    box on
    

subplot(4,1,2)
        hold on
        L_HEX  = 0.261;
        v_ion  = sqrt(2*(4200-V_bias)*1.602e-19/(18*AMU));
        T_0HEX = L_HEX/v_ion*1e9;
        tau_F2 = (F2_f_TOF(I_det1_F2)*1e9-mean_TOF2)/T_0HEX;
        rho_F2 = (F2_f_R(I_det1_F2)/(L_HEX));
        rho_over_tau_F2     = rho_F2./tau_F2;
        one_over_one_tau_F2 = 1-1./tau_F2;
        
        tau_bin = (-0.5:0.005:0.5);
        rho_bin = (-0.5:0.005:0.5);
       
        b1          = hist(one_over_one_tau_F2,tau_bin);
        rhotau_hist = histogram_2D_new(one_over_one_tau_F2',rho_over_tau_F2',tau_bin,rho_bin); 
        
        plot(tau_bin,b1,'-','LineWidth',2)
        plot(tau_bin,smooth(b1,5),'-r','LineWidth',2)
         
        xlim([-0.5 0.3]);
        %ylim([-50 50]+y_center)
        xlabel('1-1/\tau_1','FontSize',14)
        ylabel('Intensity','FontSize',14)

        set(gca,'FontSize',14)
    hold off
    box on

    
    
    subplot(4,1,3)
        hold on
      
        pcolor(tau_bin,rho_bin,rhotau_hist)
        colorbar; shading flat; colormap('hot')%default
%         R=40;t=(0:0.1:2*pi);x=R*cos(t);y=R*sin(t);plot(x,y,'-w')
         xlim([-0.2 0.2])
         ylim([ 0.0 0.2])
%         xlabel('Horizontal (mm)','FontSize',14)
%         ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
    hold off    
    
    
   
subplot(4,1,4)
        hold on
        z_center = 0;
        y_center = 0
        
        z_bin = (-50:1:50);
        y_bin = (-50:1:50) + y_center;
        %xy_hist = histogram_2D_new(F1_f_Z(I_det1_F1)*(1e3),F1_f_Y(I_det1_F1)*(1e3),z_bin,y_bin);
        xy_hist = histogram_2D_new(F2_f_Z(I_det1_F2)'*(1e3),F2_f_Y(I_det1_F2)'*(1e3),z_bin,y_bin);
        
        pcolor(z_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        R=40;t=(0:0.1:2*pi);x=R*cos(t);y=R*sin(t);plot(x,y,'-w')
        xlim([-50 50])
        ylim([-50 50]+y_center)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
         x=(-40:1:40);y=20+x*tan(60/180*pi);plot(x,y,'c') %y<20+x*alfa
         x=(-40:1:40);y=20-x*tan(60/180*pi);plot(x,y,'c') %y<20-x*alfa 
         
         x=(-40:1:40);y=-20+x*tan(60/180*pi);plot(x,y)
         x=(-40:1:40);y=-20-x*tan(60/180*pi);plot(x,y)
         
         x=(-40:1:40);y=10+x*0;plot(x,y)
         x=(-40:1:40);y=-10+x*0;plot(x,y,'c')  % y>10-x0
         
         t=(0:0.1:2*pi+0.1);x=10*cos(t);y=10*sin(t);plot(x,y)
    hold off    
    
   
    
    
end
%% %% plot of TOF for DET 2
%do_DLD=0
if (do_DLD==1)
det2_Ypos_min=-0.40;
det2_Ypos_max= 0.40;

det2_Zpos_min=-0.40;
det2_Zpos_max= 0.40;

det2_Xpos_min=0.966;
det2_Xpos_max=0.968;


 
cut_det2_F1=   (F1_f_Y>=det2_Ypos_min).*(F1_f_Y<=det2_Ypos_max)...
             .*(F1_f_Z>=det2_Zpos_min).*(F1_f_Z<=det2_Zpos_max)...
             .*(F1_f_X>=det2_Xpos_min).*(F1_f_X<=det2_Xpos_max);


cut_det2_F2=  (F2_f_Y>=det2_Ypos_min).*(F2_f_Y<=det2_Ypos_max)...
            .*(F2_f_Z>=det2_Zpos_min).*(F2_f_Z<=det2_Zpos_max)...
            .*(F2_f_X>=det2_Xpos_min).*(F2_f_X<=det2_Xpos_max);
 

      
%Loss_F1=sum(cut_det1_F1==0)/(sum(cut_det1_F1==1)+sum(cut_det1_F1==0));
%Loss_F2=sum(cut_det1_F2==0)/(sum(cut_det1_F2==1)+sum(cut_det1_F2==0));


I_det2_F1=find(cut_det2_F1);
I_det2_F2=find(cut_det2_F2);



 mean_TOF1=mean(F1_i_TOF)*1e9
 mean_TOF2=mean(F2_i_TOF)*1e9

 

 
 
figure(5);clf
    subplot(2,1,1)
        hold on
        x=(0:1:12000)
        a1=hist(F1_f_TOF(I_det2_F1)*1e9-mean_TOF1,x);
        a2=hist(F2_f_TOF(I_det2_F2)*1e9-mean_TOF2,x);
  
        plot(x,a1,'r','LineWidth',2)
        plot(x,a2,'b','LineWidth',2)
        %plot(x,a2,'b','LineWidth',3)
        xlim([1000 8000])
        xlabel('TOF (ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        %legend('0.87 eV')
        %legend('1 eV','2.3 eV')
    hold off
    box on
    

    
 subplot(2,1,2)
        hold on
        z_center = 0;
        y_center = 0
        
        z_bin   = (-50:1:50);
        y_bin   = (-50:1:50)+y_center;
        xy_hist = histogram_2D_new(F1_f_Z(I_det2_F1)'*(1e3),F1_f_Y(I_det2_F1)'*(1e3),z_bin,y_bin);
        %xy_hist=histogram_2D_new(F2_f_Z(I_det2_F1)'*(1e3),F2_f_Y(I_det2_F1)'*(1e3),z_bin,y_bin);
        
        pcolor(z_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        R=40;t=(0:0.1:2*pi);x=R*cos(t);y=R*sin(t);plot(x,y,'-w')
        xlim([-50 50])
        ylim([-50 50]+y_center)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
%          x=(-40:1:40);y=20+x*tan(60/180*pi);plot(x,y,'c') %y<20+x*alfa
%          x=(-40:1:40);y=20-x*tan(60/180*pi);plot(x,y,'c') %y<20-x*alfa 
%          
%          x=(-40:1:40);y=-20+x*tan(60/180*pi);plot(x,y)
%          x=(-40:1:40);y=-20-x*tan(60/180*pi);plot(x,y)
%          
%          x=(-40:1:40);y=10+x*0;plot(x,y)
%          x=(-40:1:40);y=-10+x*0;plot(x,y,'c')  % y>10-x0
         
         %t=(0:0.1:2*pi+0.1);x=10*cos(t);y=10*sin(t);plot(x,y)
    hold off    
    
  
   
    
    
end

 
%% %% plot of TOF for DET 3
%do_HCH=1
if (do_HCH==1)

%%% geometrical cuts definition    
    det3_Ypos_min=0.125;
    det3_Ypos_max=0.195;

    det3_Zpos_max=0.040;
    det3_Zpos_min=-0.040;

    det3_Xpos_min=0.583;
    det3_Xpos_max=0.623;


 
cut_det3_F1=  (F1_f_Y>=det3_Ypos_min).*(F1_f_Y<=det3_Ypos_max)...
             .*(F1_f_Z>=det3_Zpos_min).*(F1_f_Z<=det3_Zpos_max)...
             .*(F1_f_X>=det3_Xpos_min).*(F1_f_X<=det3_Xpos_max);

cut_det3_F2=  (F2_f_Y>=det3_Ypos_min).*(F2_f_Y<=det3_Ypos_max)...
          .*(F2_f_Z>=det3_Zpos_min).*(F2_f_Z<=det3_Zpos_max)...
          .*(F2_f_X>=det3_Xpos_min).*(F2_f_X<=det3_Xpos_max);
 

cut_det1_F1= (F1_f_X>=0.262).*(F1_f_X<=0.264);
cut_det1_F2= (F2_f_X>=0.262).*(F2_f_X<=0.264);

      
Loss_F1=sum(cut_det3_F1==0)/(sum(cut_det3_F1==1)+sum(cut_det3_F1==0));
Loss_F2=sum(cut_det3_F2==0)/(sum(cut_det3_F2==1)+sum(cut_det3_F2==0));

Loss_F2_hEX=sum(cut_det1_F2==1)/(sum(cut_det3_F2==1)+sum(cut_det3_F2==0));

clear E_k Trans_E_F1 Trans_E_F2
for jj=1:10
Ejj=jj;
cut_det3_E_F1=(F1_E_release==Ejj);
cut_det3_E_F2=(F2_E_release==Ejj);

    cut_det3_all_F1=cut_det3_F1.*cut_det3_E_F1;
    cut_det3_all_F2=cut_det3_F2.*cut_det3_E_F2;
    E_k(jj)=Ejj;
    Trans_E_F1(jj)=sum(cut_det3_all_F1==1)/(sum(cut_det3_E_F1==1));
    Trans_E_F2(jj)=sum(cut_det3_all_F2==1)/(sum(cut_det3_E_F2==1));
end

cut_det3_E_F1=(F1_E_release>0).*(F1_E_release<40);
cut_det3_E_F2=(F2_E_release>0).*(F2_E_release<40);

I_det3_F1=find(cut_det3_F1.*cut_det3_E_F1);
I_det3_F2=find(cut_det3_F2.*cut_det3_E_F2);

text_1=strcat('M1, T =',num2str(1-Loss_F1))
text_2=strcat('M2, T =',num2str(1-Loss_F2))

 mean_TOF1=mean(F2_i_TOF)*1e9;
 mean_TOF2=mean(F2_i_TOF)*1e9;

 
figure(2);clf
title('Detection on DET 3')
hold on
    h1=plot(E_k,Trans_E_F1,'.r')
    h2=plot(E_k,Trans_E_F2,'.b')
    ylim([0 1.1])
     xlabel('E_k (eV)','FontSize',14)
     ylabel('Tranmission','FontSize',14)
     set(gca,'FontSize',14)
    legend([h1 h2],'Fragment 1','Fragment 2','FontSize',14)
    legend boxoff
hold off
box on
 
 
figure(3);clf
    subplot(3,1,1)
        hold on
        x=(0:2:12000);
        a1=hist(F1_f_TOF(I_det3_F1)*1e9-mean_TOF2,x);
        a2=hist(F2_f_TOF(I_det3_F2)*1e9-mean_TOF2,x);
  
        plot(x,a1,'r','LineWidth',2)
        plot(x,a2,'b','LineWidth',2)
        plot(x,a1+a2,'--k','LineWidth',1)
        
        %plot(x,a2,'b','LineWidth',3)
        xlim([1000 7500]+000)
        xlabel('TOF (ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        legend(text_1,text_2)
        %legend('1 eV','2.3 eV')
    hold off
    box on

subplot(3,1,2)
        hold on
        x_center=603
        z_center=0;
        y_center=160% 140
        z_bin=(-100:1:50);
        y_bin=(-100:1:50)+y_center*0;
        theta = 60/180*pi;
        X_det3_F1 = cos(theta)*(F1_f_X*1e3-x_center) - sin(theta)*(F1_f_Y*1e3-y_center);
        Y_det3_F1 = sin(theta)*(F1_f_X*1e3-x_center) + cos(theta)*(F1_f_Y*1e3-y_center);
        Z_det3_F1 = F1_f_Z*1e3;
        
        xy_hist=histogram_2D_new(Z_det3_F1(I_det3_F1)',Y_det3_F1(I_det3_F1)',z_bin,y_bin);
        
        
        %xy_hist=histogram_2D_new(F1_f_Z(I_det3_F1)'*(1e3),F1_f_Y(I_det3_F1)'*(1e3),z_bin,y_bin);
        %%xy_hist=histogram_2D_new(F2_f_Z(I_det3_F2)*(1e3),F2_f_Y(I_det3_F2)*(1e3),z_bin,y_bin);
        
        pcolor(z_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        R=40;t=(0:0.1:2*pi);x=R*cos(t);y=R*cos(30/180*pi)*sin(t)+y_center*0;plot(x,y,'-r')
        xlim([-50 50])
        ylim([-50 50]+y_center*0)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)
        legend('M1')
        set(gca,'FontSize',14)
    hold off
    box on

    
    subplot(3,1,3)
        hold on
        x_center=603
        z_center=0;
        y_center=160% 140
        z_bin=(-50:1:50);
        y_bin=(-50:1:50)+y_center*0;
        
        X_det3_F2 = cos(theta)*(F2_f_X*1e3-x_center) - sin(theta)*(F2_f_Y*1e3-y_center);
        Y_det3_F2 = sin(theta)*(F2_f_X*1e3-x_center) + cos(theta)*(F2_f_Y*1e3-y_center);
        Z_det3_F2 = F2_f_Z*1e3;
        
        xy_hist=histogram_2D_new(Z_det3_F2(I_det3_F2)',Y_det3_F2(I_det3_F2)',z_bin,y_bin);
        
        
        %xy_hist=histogram_2D_new(F1_f_Z(I_det3_F1)*(1e3),F1_f_Y(I_det3_F1)*(1e3),z_bin,y_bin);
        %xy_hist=histogram_2D_new(F2_f_Z(I_det3_F2)'*(1e3),F2_f_Y(I_det3_F2)'*(1e3),z_bin,y_bin);
        
        pcolor(z_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        R=40;t=(0:0.1:2*pi);x=R*cos(t);y=R*cos(30/180*pi)*sin(t)+y_center*0;plot(x,y,'-r')
        xlim([-50 50])
        ylim([-50 50]+y_center*0)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)
        legend('M2')
        set(gca,'FontSize',14)
    hold off
    box on
    
    
     mean(F1_f_TOF*1e9-mean_TOF1)
     mean(F2_f_TOF*1e9-mean_TOF2)
 
    
    To_DET3 = L_DET3/mean(Ion_Vtotal);
    
    F1_f_one_minus_one_over_tau = 1-To_DET3./(F1_f_TOF-mean_TOF1*1e-9);
    F2_f_one_minus_one_over_tau = 1-To_DET3./(F2_f_TOF-mean_TOF2*1e-9);
    
    mean(F1_f_one_minus_one_over_tau)
    mean(F2_f_one_minus_one_over_tau)
     
    figure(6);clf
    hold on
        x=(-1:0.0001:1);
        a1=hist(F1_f_one_minus_one_over_tau(I_det3_F1),x);
        a2=hist(F2_f_one_minus_one_over_tau(I_det3_F2),x);
  
        plot(x,a1,'r','LineWidth',2)
        plot(x,a2,'b','LineWidth',2)
        %plot(x,a2,'b','LineWidth',3)
        %xlim([1000 8000])
        xlabel('TOF (ns)','FontSize',14)
    
    hold off
    
    
end   
    
    
    

%%%%ELECTRON SPECTROMTER ANALYSIS START %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% plot of TOF and XY for det-UP
%do_ELEC=0
if (do_ELEC==1)
E_detup_Ypos=0.142  ;% 0.291;     %-0.067 ;       %0.212; 
E_detdw_Ypos=-0.142;%0.009%    %-0.067 ;       %0.212;
 
 cut_det_up=(F2_f_Y>=E_detup_Ypos*0.99).*(F2_f_Y<=E_detup_Ypos*1.01).*(F2_E_release<40);
 cut_det_dw=(F2_f_Y<=E_detdw_Ypos*0.99).*(F2_f_Y>=E_detdw_Ypos*1.01).*(F2_E_release<40);
 
 %cut_det_dw=(F2_f_Y>=E_detdw_Ypos*0.99).*(F2_f_Y<=E_detdw_Ypos*1.01).*(F2_E_release<40); %%%NB
 
 
 cut_det_all=cut_det_up+cut_det_dw;
 
 
 E_vec = [0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 4.3 4.8];

 E1=20;     E2=0.8;   E3=1.8;    E4=4.2;      E5=4.2;
 
 %%% UP
 I_detup=find(cut_det_up==1);
 I_detup_lost=find(cut_det_up==0);
 I_detup_E1=find((cut_det_up==1).*(F2_E_release==E1));
 I_detup_E2=find((cut_det_up==1).*(F2_E_release==E2));
 I_detup_E3=find((cut_det_up==1).*(F2_E_release==E3));
 I_detup_E4=find((cut_det_up==1).*(F2_E_release==E4));
 I_detup_E5=find((cut_det_up==1).*(F2_E_release==E5));
 %%%%%
% 
%  for jj_E=1:[2 4 5 7 9]
%         E_now=E_vec(jj_E);
%         eval(strcat('I_detup_E',num2str(jj_E),'=find((cut_det_up==1).*(F2_E_release==',num2str(E_now),'));'));
%  end
 %%%% DW
 %E1=6;      E2=4;       E3=2;
%  I_detdw=find(cut_det_dw==1);
%  I_detdw_lost=find(cut_det_dw==0);
%  I_detdw_E1=find((cut_det_dw==1).*(F2_E_release==E1));
%  I_detdw_E2=find((cut_det_dw==1).*(F2_E_release==E2));
%  I_detdw_E3=find((cut_det_dw==1).*(F2_E_release==E3));
 %%%%%%
 
 
 
 lost        = sum(cut_det_all==0)/( sum(cut_det_all==0) + sum(cut_det_all==1) )  ;
 transmision = sum(cut_det_all==1)/( sum(cut_det_all==0) + sum(cut_det_all==1) )  ; 

mean_TOF2=mean(F2_i_TOF)*1e9; % not the correct time zero for the transversely emitted electrons


size(I_detup)

figure(2);clf

    %subplot(2,1,1)
    text=strcat('transmission =',num2str(transmision))
    title(text)
        hold on
        x=(0:2:1000);
  w=0       
            d=F2_f_TOF(I_detup_E1)*1e9-F2_i_TOF(I_detup_E1)*1e9;
                noise=(rand(1,length(d),'double')-0.5)*w;
                a1=hist(d+noise',x);
            h1=stairs(x,a1,'Color','b','LineWidth',1.5);
  
            d=F2_f_TOF(I_detup_E2)*1e9-F2_i_TOF(I_detup_E2)*1e9;
                noise=(rand(1,length(d),'double')-0.5)*w;
                a1=hist(d+noise',x);
            h2=stairs(x,a1,'Color','r','LineWidth',1.5);
  
            d=F2_f_TOF(I_detup_E3)*1e9-F2_i_TOF(I_detup_E3)*1e9;
                noise=(rand(1,length(d),'double')-0.5)*w;
                a1=hist(d+noise',x);
            h3=stairs(x,a1,'Color','k','LineWidth',1.5);
  
            d=F2_f_TOF(I_detup_E4)*1e9-F2_i_TOF(I_detup_E4)*1e9;
                noise=(rand(1,length(d),'double')-0.5)*w;
                a1=hist(d+noise',x);
            h4=stairs(x,a1,'Color','g','LineWidth',1.5);
 
            d=F2_f_TOF(I_detup_E5)*1e9-F2_i_TOF(I_detup_E5)*1e9;
                noise=(rand(1,length(d),'double')-0.5)*w;
                a1=hist(d+noise',x);
            h5=stairs(x,a1,'Color','c','LineWidth',1.5);
%  
%             
%   
%         for jj_plot=1:10
%             jj_plot
%             eval(strcat('I_now=I_detup_E',num2str(jj_plot),';'));
%             E_plot=E_vec(jj_plot);
%             d=F2_f_TOF(I_now)*1e9-F2_i_TOF(I_now)*1e9;
%             noise=(rand(1,length(d),'double')-0.5)*w;
%             a1=hist(d+noise',x);
% 
%             h0=stairs(x,a1,'Color',[0.1 0.1 0.1],'LineWidth',1.5);
%             eval(strcat('h',num2str(jj_plot),'=h0'));
%             txt1=strcat(num2str(E_plot),blanks(1),' eV')
%             eval(strcat('legend([h',num2str(jj_plot),'],txt1)'))
%         end
%         
       
        
        xlim([40 2000])
        xlabel('TOF Det UP (ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        %legend('0.87 eV')
        
        
         txt1=strcat(num2str(E1),blanks(1),' eV');
         txt2=strcat(num2str(E2),blanks(1),' eV');
         txt3=strcat(num2str(E3),blanks(1),' eV');
         txt4=strcat(num2str(E4),blanks(1),' eV');
         txt5=strcat(num2str(E5),blanks(1),' eV');
         legend([h1 h2 h3 h4 h5],txt1,txt2,txt3,txt4,txt5)
    hold off
    box on

figure(22);clf
%subplot(2,1,2)
        hold on
        x_center=0
        x_bin=(-40:1:40);
        y_bin=(-40:1:40)+x_center;
        %I_ana=I_detup;
        I_ana=I_detup_E1;
        
        xy_hist=histogram_2D_new(F2_f_X(I_ana)'*(1e3),F2_f_Z(I_ana)'*(1e3),x_bin,y_bin);
        
        pcolor(x_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        R=20;t=(0:0.1:2*pi);x=R*cos(t);y=R*sin(t)+x_center;plot(x,y,'-r')
        xlim([-20 20])
        ylim([-20 20]+x_center)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
    hold off
    box on

end
    return
%% save splines to mat file
% save('simulation.mat','PP_28eV','PP_24eV','PP_15eV')
% figure(500);clf
% hold on
%     %a1=hist(F2_f_TOF(I_detup_E1)*1e9-mean_TOF2,x);
%      x=(0:2.0:1000);
%     a1=hist(F2_f_TOF(I_detup_E1)*1e9-F2_i_TOF(I_detup_E1)*1e9,x);
%     stairs(x,a1,'--b','LineWidth',1.5)
%     xlim([0 200])
%     
%     PP_28eV=spline(x,a1)
%     a_sp=ppval(PP_28eV,x)
%  
%     plot(x,a_sp,'r')
% hold off    
%%%%ELECTRON SPECTROMTER ANALYSIS STOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 




return 
%% plot of TOF and XY for det 1 (HEX detector)

 E_det1_Xpos= 0.26100; %0.217 ;  %0.212;
 
 cut_det2=(F2_f_X>E_det1_Xpos-0.01).*(F2_f_X<E_det1_Xpos+0.01);
 I_det2=find(cut_det2);
 I_det2_lost=find(cut_det2==0);

 
 

mean_TOF1=mean(F2_i_TOF)*1e9;
size(I_det2)
figure(2);clf
    subplot(2,1,1)
        hold on
        x=(0:10:2500)
        a1=hist(F2_f_TOF(I_det2)*1e9-mean_TOF1,x);
       
        plot(x,a1,'r','LineWidth',3)
        %plot(x,a2,'b','LineWidth',3)
        xlim([00 2500])
        xlabel('TOF after laser(ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        %legend('0.87 eV')
        %legend('1 eV','2.3 eV')
    hold off
    box on
subplot(2,1,2)
        hold on

        x_bin=(-40:1:40);
        y_bin=(-40:1:40);
        xy_hist=histogram_2D_new(F2_f_Z(I_det2)*(1e3),F2_f_Y(I_det2)*(1e3),x_bin,y_bin);
        
        pcolor(x_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        t=(0:0.1:2*pi);x=40*cos(t);y=40*sin(t);plot(x,y,'-r')
        xlim([-40 40])
        ylim([-40 40])
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
    hold off
    box on

%% stop
return
  
%% plot of TOF and XY for det 3

E_det3_Xpos=0.465; %-0.067 ;  %0.212; 
E_det3_Ypos=0.0; %-0.067 ;  %0.212;
 
 cut_det3=(F1_f_X>E_det3_Xpos).*(F2_f_Y<E_det3_Ypos);
 I_det3=find(cut_det3);
 I_det3_lost=find(cut_det3==0);

 
 

mean_TOF2=mean(F1_i_TOF)*1e9;
size(I_det2)
figure(2);clf
    subplot(2,1,1)
        hold on
        x=(0:10:10000)
        a1=hist(F1_f_TOF(I_det3)*1e9-mean_TOF2,x);
       
        plot(x,a1,'r','LineWidth',2)
        %plot(x,a2,'b','LineWidth',3)
        xlim([5500 6500])
        xlabel('TOF after laser(ns)','FontSize',14)
        %set(gca,'FontSize',14,'Xtick',[0 50 100 150 200 250 300 350 400  ])
        %legend('0.87 eV')
        %legend('1 eV','2.3 eV')
    hold off
    box on

subplot(2,1,2)
        hold on
        x_center=-100
        x_bin=(-40:1:40);
        y_bin=(-40:1:40)+x_center;
        xy_hist=histogram_2D_new(F1_f_Z(I_det3)*(1e3),F1_f_Y(I_det3)*(1e3),x_bin,y_bin);
        
        pcolor(x_bin,y_bin,xy_hist)
        colorbar; shading flat; colormap('hot')%default
        t=(0:0.1:2*pi);x=40*cos(t);y=40*sin(t)+x_center;plot(x,y,'-r')
        xlim([-40 40])
        ylim([-40 40]+x_center)
        xlabel('Horizontal (mm)','FontSize',14)
        ylabel('Vertical (mm)','FontSize',14)

        set(gca,'FontSize',14)
    hold off
    box on



    
    
return    
%%

return
%%
 figure(2);clf
    hold on
    plot(Ion_R*(1e3),Ion_Vr/mm_microsec2m_sec)
 
 
 
   subplot(2,2,1)
        hist(F1_i_Vx,1000)
   
   subplot(2,2,2)
        hist(F1_i_Y,(-5:0.3:5))
   
   subplot(2,2,3)
        hist(F1_i_Z,(-5:0.3:5))
   
   
 
  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear X_Edet_bin Z_Edet_bin XZ_hist_ion
 X_Edet_bin=(0:80)*0.5-20;
 Y_Edet_bin=(0:80)*0.5-20;
 Z_Edet_bin=(0:80)*0.5-20;
 
    XZ_hist_Edet=histogram_2D_new(F2_f_X(I_2_Dup)*(1e3),F2_f_Z(I_2_Dup)*(1e3),X_Edet_bin,Z_Edet_bin);   
 
    
   
    
    
    figure(10);clf
        colormap('default');
        hold on
        pcolor(X_Edet_bin,Z_Edet_bin,XZ_hist_Edet);
        colorbar
        shading flat
        ylabel('Ions horizontal coordinate Z(mm)','FontSize',10)
        xlabel('Ions longitudinal X(mm)','FontSize',10)
        hold off
        box on
 
print -depsc2 figure10.ps
 
   
 return       
        
%%        
    figure(11);clf
        hold on
        hist((F2_f_R),100)
        
        
   
    clear R_Edet_bin TOF_Edet_bin RT_hist_Edet
     R_Edet_bin  =(0:80)*0.25;
     TOF_Edet_bin=(0:250)*0.25;
     RT_hist_Edet=histogram_2D_new(F2_f_TOF(I_2_Dup)*(1e9),F2_f_R(I_2_Dup)*(1e3),TOF_Edet_bin,R_Edet_bin);
    figure(12);clf
        colormap('default');
        hold on
        pcolor(TOF_Edet_bin,R_Edet_bin,RT_hist_Edet);
        colorbar
        shading flat
        axis([0 60 0 20])
        ylabel('R(mm)','FontSize',16)
        xlabel('TOF(ns)- TOF_{offset}','FontSize',16)
        set(gca,'FontSize',16)
        hold off
        box on
    
    print_pdf('','figure12')
        
        
    figure(13);clf
        hold on
        theta_bin=(0:0.02:pi)
        a=hist(F2_i_theta(I_2_lost),theta_bin)
        b=hist(F2_i_theta(I_2_Dup),theta_bin)
        
        plot(theta_bin,a,'-b')
        plot(theta_bin,b,'-r')
        %plot(theta_bin,a-b,'-g')
        hold off
    
        
   figure(14);clf
        subplot(2,1,1)
            hold on
            %hist(F2_i_phi_xz(I_2_Dup),(0:0.1:2*pi))
            plot(F2_i_phi_xz(I_2_Dup),F2_f_R(I_2_Dup)*1000,'.')
            %plot(F2_i_phi,F2_f_phi,'.')
             xlabel('\phi XZ initial')
            ylabel('R final')
            hold off
        subplot(2,1,2)
            plot(F2_i_phi_xz(I_2_Dup),F2_f_phi_xz(I_2_Dup),'.')
            xlabel('\phi XZ initial')
            ylabel('\phi XZ final')
            
            
   figure(15);clf
        hold on
         plot(F2_i_phi_xz(I_2_Dup),F2_i_theta(I_2_Dup)*1000,'.')
         %plot(F2_f_phi_xz(I_2_Dup),F2_f_R(I_2_Dup)./F2_f_TOF(I_2_Dup)*1000,'.')
        %axis([0 2*pi 0 1])
        hold off
        
    figure(16);clf
        hold on
        plot(F2_i_phi_xz(I_2_Dup),F2_f_TOF(I_2_Dup),'.')
        
        hold off
    
        
        
        
        
    figure(17);clf
        subplot(2,1,1)
            hold on
            plot(F2_f_TOF(I_2_Dup)*1e9,F1_f_TOF(I_2_Dup),'.')
            hold off
            
        subplot(2,1,2)
            hold on
            plot(F2_f_R(I_2_Dup)*1e9,F1_f_R(I_2_Dup),'.')
            hold off
            
            