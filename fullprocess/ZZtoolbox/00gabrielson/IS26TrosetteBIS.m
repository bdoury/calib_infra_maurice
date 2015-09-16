function transf = IS26TrosetteBIS(freq_vector,NRS)
%==========================================================================
% Thermo-viscous acoustical model for a rosette pipe system
%
% This version contains an example system based on the L4
% element, one of the 70-m rosette elements at I57US. It is easily
% modified for other rosette systems.
%
% T. Gabrielson (tbg3@psu.edu)
% Penn State University
% [January 2012]
%
%==========================================================================
%
%=============== USER INPUT BLOCK =======================================
% Pipe specs
%  Pipe specs:
%        numbers and lengths [m] of primary, secondary, and tertiary pipes

Npri  = NRS.Npri;      
Nsec  = NRS.Nsec;    

Lpri  = NRS.Lpri;  
Lsec  = NRS.Lsec;  

a_pri = NRS.a_pri  ; %  DN12
a_sec = NRS.a_sec ;  %  DN12

a_sec_cavity = NRS.a_sec_cavity;
L_sec_cavity = NRS.L_sec_cavity;

a_microbarometer = NRS.a_microbarometer;  %  radius and length to get correct MB2000
L_microbarometer = NRS.L_microbarometer;  %  volume (0.0006 m^3 from Alcoverro)

L_stub = NRS.L_stub;
a_stub = NRS.a_stub;
%
%
if L_stub == 0
    a_pri_cavity = a_microbarometer;
    L_pri_cavity = L_microbarometer;
else
    a_pri_cavity = a_sec_cavity;
    L_pri_cavity = L_sec_cavity;
end

%==============================================
RS_at_primary_manifold   = false;
RS_at_secondary_manifold = false;
%   Resonance suppressor specs
%        If either of the 'RS_at_..._manifold' flags are false, the
%        corresponding radius and length inputs are ignored.
a_rsp = 0.975e-3;  
L_rsp = 20e-3;  %  Radius and length of capillary [m]
a_rss = 0.975e-3;  
L_rss = 20e-3;  %  Radius and length of capillary [m]
%   RS_xxx means that the suppressors are at the
%       outboard ends of the primary pipes.
%   RS_at_primary_manifold means that a suppressor is at the
%       primary-manifold end of the stub to the microbarometer.
%==============================================
%   Temperature [deg C] and ambient pressure [Pa]
tempC = 5;  
atmos_press = 88000;  %   88000 Pa is roughly "standard"
% atmospheric pressure at the altitude of
% Conrad Observatory
%==============================================
% 'damping_correction' is an empirical correction factor used in
%  Tmatrix_air.  Set to zero to disable that correction.
damping_correction = 0;  %  normally set to one (???)
%==============================================
%  Miscellaneous pre-calculations
tempK     = tempC + 273;  
cc        = 20.05*sqrt(tempK); % sound celerity 
MW_air    = 28.97;  
Rgas      = 8314.3;  
rho       = atmos_press*MW_air/(Rgas*tempK);
%
% The ultimate model result is the ratio of (complex) pressure at the
%   microbarometer to acoustic pressure at the inlet:  p_ratio
transf = zeros(size(freq_vector)); 
%
% Set length increment for pipe model based on highest frequency
% i.e.the smallest wavelength
freq_max       = max(freq_vector);
min_wavelength = cc/freq_max;
max_dx         = min_wavelength/1000;
%
% Set some transfer matrices to the identity matrix.  This permits
% omitting the feature (by leaving the matrix as the identity matrix or
% including the feature by changing the values.
TT_RS_sec  = [1 0; 0 1];   %  resonance suppressor at sec manifold
TT_RS_pri  = [1 0; 0 1];   %  resonance suppressor at pri manifold
TT_stub    = [1 0; 0 1];   %  stub to microbarometer
%
%   Primary computational loop (one pass per frequency)
%
for ii = 1:length(freq_vector)
    freq = freq_vector(ii);
    %
    %  Transfer matrices for pipe sections
    %
    %         Primary pipe
    Nslice_pri = floor(Lpri/max_dx) + 1;  %  number of "slices"
    TT_pri_item = Tmatrix_air(freq, a_pri, Lpri/Nslice_pri, tempK,...
        atmos_press, damping_correction); % transfer matrix for pipe slice
    TT_pri = TT_pri_item^Nslice_pri;  %  transfer matrix for entire length of pipe
    %
    %         Secondary pipe
    Nslice_sec = floor(Lsec/max_dx) + 1;
    TT_sec_item = Tmatrix_air(freq, a_sec, Lsec/Nslice_sec, tempK,...
        atmos_press, damping_correction);
    TT_sec = TT_sec_item^Nslice_sec;
    %
    %         Central-manifold-to-microbarometer pipe (if used)
    if L_stub > 0  %  otherwise, leave matrix as identity matrix
        Nslice_stub = floor(L_stub/max_dx) + 1;
        TT_stub_item = Tmatrix_air(freq, a_stub, L_stub/Nslice_stub, ...
            tempK, atmos_press, damping_correction);
        TT_stub = TT_stub_item^Nslice_stub ;
    end
    %
    %  Radiation admittance and corresponding transfer matrix
    Y_rad  = radiation_admittance(rho, cc, a_sec, freq);
    TT_rad = [1  1/Y_rad;  0  1];
    %
    Y_cavity_sec     = Y_cavity(freq, a_sec_cavity, L_sec_cavity,...
        tempK, atmos_press);
    Y_cavity_pri     = Y_cavity(freq, a_pri_cavity, L_pri_cavity,...
        tempK, atmos_press);
    %
    Y_microbarometer = Y_cavity(freq, a_microbarometer,...
        L_microbarometer, tempK, atmos_press);
    %
    %  Resonance suppressors if used
%     if RS_at_secondary_manifold  %  otherwise, matrix is identity matrix
%         Zrs = Z_capillary_inlet(freq, a_rss, L_rss, tempK, atmos_press);
%         TT_RS_sec = [1 Zrs; 0 1];
%     end
%     %
%     if RS_at_primary_manifold  %  otherwise, matrix is identity matrix
%         Zrs = Z_capillary_inlet(freq, a_rsp, L_rsp, tempK, atmos_press);
%         TT_RS_pri = [1 Zrs; 0 1];
%     end
    %
    %     Total "passive" (i.e., all inlet pressures = 0)
    %               tertiary-manifold admittance
    %
    %     Translate passive tertiary admittance to secondary manifold
    Y_rad_xlate = admittance_transfer(TT_sec, Y_rad);
    Y_total_sec = Y_cavity_sec + Nsec*Y_rad_xlate;
    %
    %     Translate passive secondary admittance to primary manifold
    Y_sec_xlate = admittance_transfer(TT_pri, Y_total_sec);
    Y_total_pri = Y_cavity_pri + (Npri - 1)*Y_sec_xlate;
    %
    %     Translate central (primary-total) admittance outbound along
    %               primary pipe
    Y_outpri_xlate = admittance_transfer(TT_RS_sec*TT_pri, Y_total_pri);
    %
    %     Add Nsec - 1 secondary pipes and secondary cavity
    Y_sec_reference = Y_outpri_xlate + Y_cavity_sec +...
        (Nsec - 1)*Y_rad_xlate;
    %
    %
    %
        %
%     % Compute pressure ratio, p_in/p_sec
%     p_rat_sec = [1 0]*TT_sec*[1; Y_sec_reference];
%     %
%     % Compute pressure ratio, p_sec/p_pri
%     p_rat_pri = [1 0]*(TT_pri_OB*TT_suppr*TT_pri_IB)*[1; Y_total_pri];
%     %
%     % Extra feed pipe to microbarometer: pressure ratio
%     p_rat_stub = [1 0]*TT_stub*[1; Y_microbarometer];
%     %
%     % Compute the pressure ratio, p_pri/p_in
%     transf(ii) = 1/(p_rat_sec*p_rat_pri*p_rat_stub);
%     %
% 
    %================================================
    %
    %  Compute pressure ratio, p_inlet/p_ter
    %
    %  Compute pressure ratio, p_ter/p_sec
    p_ratio_sec  = [1 0]*TT_rad*TT_sec*[1; Y_sec_reference];
    %
    %  Compute pressure ratio, p_sec/p_pri
    p_ratio_pri  = [1 0]*TT_pri*[1; Y_total_pri];
    %
    %  Extra feed pipe to microbarometer: pressure ratio
    p_ratio_stub = [1 0]*TT_stub*[1; Y_microbarometer];
    %
    %  Compute the pressure ratio, p_pri/p_in
    transf(ii) = 1/(p_ratio_sec*p_ratio_pri*p_ratio_stub);
    %
end
    %
end
%   Sum over all inlets (equivalent to vertical incidence)
% p_ratio = p_ratio;
