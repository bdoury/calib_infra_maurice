function [p_totalNRSsensor, freq_vector, p_total_NRS, TF_ref_sensor, ...
    TF_ref_sensor4freqRatio] = ...
    theoreticalSUTresponse4IS26(...
    allfrqsPfilters, sensor_UT, ref_sensor,firflag)
%
%
%
% used functions: idc2fap_is.m (from Benoit Doury)
%                 PW_70_meter_I57.m

%==========================================================================
% Thermo-viscous acoustical model for a rosette pipe system
%
% T. Gabrielson (tbg3@psu.edu)
% Penn State University
% [January 2012]
%
%==========================================================================
addpath ZZtoolbox/00gabrielson/
addpath ZZtoolbox/00benoit/

% elevation has no effect on small antenna as IS26
% Indeed the min wavelength is around 300/5 = 60 m much larger
% than the SUT aperture.
arrival_theta_degrees = 0; % Elevation angle: horizontal = 0
arrival_theta_radians = arrival_theta_degrees*pi/180;
%=============== USER INPUT BLOCK =======================================
% Description of IS26
%
% Pipe specs
% Npri = 8;  % number of primary pipes
% Nsec = 18; % number of secondary pipes
% Lpri = 27; % lengths [m] of primary pipes
% Lsec = 8;  % lengths [m] of secondary pipes
Npri   = 4;  % number of primary pipes
Nsec   = 24; % number of secondary pipes
Lpri   = 5.75; % lengths [m] of primary pipes
Lsec   = 3.25;  % lengths [m] of secondary pipes
%
% Primary and secondary pipe radii [m]
a_pri  = 0.824*25.4e-3/2; % Sch40 galvanized iron pipe (3/4 inch)
a_sec  = 0.622*25.4e-3/2; % Sch40 galvanized iron pipe (1/2 inch)
%
% Length of stub from primary manifold to microbarometer [m]
Lstub  = 1.6;%= 3.708;
a_stub = a_sec;
%
% Cavity specs (primary and secondary manifolds)
%
a_pri_cavity = 5/100;   % 5 cm radius
L_pri_cavity = 3.8/100; % 3.8 cm height
a_sec_cavity = a_pri_cavity;
L_sec_cavity = L_pri_cavity;
% radius and length to get correct MB2000
a_microbarometer = 7/100;
L_microbarometer = 3.9/100; % volume (0.0006 m^3 from Alcoverro)
%
% Resonance suppressor specs (ignored if resonance_suppressed = false)
resonance_suppressed = false;
a_rs = 0.95e-3;
L_rs = 20e-3; % for L3/L4 (changed below for L1/L2)
x_suppressor = Lpri;
% location of suppressor along primary pipe
% (Set to Lpri if no suppressor or if suppressor is
% in the secondary manifold.)
%
% Conditions selected for a sample measurement at I57US (Pinon Flats)
tempC    = 5;
tempK    = tempC + 273;
pressure = 88000;
rho      = 1.0;
c        = 20.05*sqrt(tempK);
%
% Select frequencies for solution
N_freq_vector = 400;
freq_vector = logspace(log10(0.001),log10(30),N_freq_vector) .';
%
% plot_f_lo = 4e-3;
% plot_f_hi = 40; % range of frequencies to display
%
% 'damping_correction' is an empirical correcction factor used in
% Tmatrix_air. Set to zero to disable that correction.
damping_correction = 0; % normally set to one
%
% Lengths of primary pipe section outboard (OB) and inboard (IB) of
% resonance suppressor
Lpri_OB = Lpri - x_suppressor;
Lpri_IB = x_suppressor;
%
% The ultimate model result is the ratio of (complex) pressure at the
% microbarometer to acoustic pressure at the inlet: p_rat
p_rat = zeros(size(freq_vector)); % Allocate space for pressure ratio
%
% Set length increment for pipe model based on highest frequency
min_wavelength = c/freq_vector(end);
max_dx         = min_wavelength/1000;
%
% Set some transfer matrices to the identity matrix. This permits
% omitting the feature (by leaving the matrix as the identity matrix or
% including the feature by changing the values.
TT_suppr = [1 0; 0 1];
% resonance suppressor
TT_pri_OB = [1 0; 0 1]; % primary pipe section outboard of suppressor
TT_stub   = [1 0; 0 1];
% stub from primary manifold to microbarometer
%
% Primary computational loop (one pass per frequency)
%
for ii = 1:N_freq_vector
    freq = freq_vector(ii);
    %
    % Transfer matrices for pipe sections
    %
    % Primary pipe inboard of resonance suppressor
    Nslice_pri = floor(Lpri_IB/max_dx) + 1; % number of "slices"
    TT_pri     = Tmatrix_air(freq, a_pri, Lpri_IB/Nslice_pri, tempK,...
        pressure, damping_correction); % transfer matrix for pipe slice
    TT_pri_IB = TT_pri^Nslice_pri; % transfer matrix for entire length of pipe
    %
    % Primary pipe outboard of resonance suppressor (if used)
    if Lpri_OB > 0 % otherwise, leave matrix as identity matrix
        Nslice_res = floor(Lpri_OB/max_dx) + 1;
        TT_res = Tmatrix_air(freq, a_pri, Lpri_OB/Nslice_res, tempK,...
            pressure, damping_correction);
        TT_pri_OB = TT_res^Nslice_res;
    end
    %
    % Secondary pipe
    Nslice_sec = floor(Lsec/max_dx) + 1;
    TT_sec1 = Tmatrix_air(freq, a_sec, Lsec/Nslice_sec, tempK,...
        pressure, damping_correction);
    TT_sec = TT_sec1^Nslice_sec;
    %
    % Central-manifold-to-stub pipe (if used)
    if Lstub > 0 % otherwise, leave matrix as identity matrix
        Nslice_stub = floor(Lstub/max_dx) + 1;
        TT_stub1 = Tmatrix_air(freq, a_stub, Lstub/Nslice_stub, tempK,...
            pressure, damping_correction);
        TT_stub = TT_stub1^Nslice_stub;
    end
    %
    % Radiation admittance
    Y_rad  = radiation_admittance(rho, c, a_sec, freq);
    TT_rad = radiation_transferM(rho, c, a_sec, freq);
    %
    % Cavity admittances (assume cylindrical cavities)
    Y_cavity_sec = Y_cavity(freq, a_sec_cavity, L_sec_cavity,...
        tempK, pressure);
    Y_cavity_pri = Y_cavity(freq, a_pri_cavity, L_pri_cavity,...
        tempK, pressure);
    Y_microbarometer = Y_cavity(freq, a_microbarometer, L_microbarometer,...
        tempK, pressure);
    %
    % Resonance suppressor if used
    if resonance_suppressed % otherwise, leave matrix as identity matrix
        Zrs = Z_capillary_inlet(freq, a_rs, L_rs, tempK, pressure);
        TT_suppr = [1 Zrs; 0 1];
    end
    %
    % Total "passive" (i.e., all inlet pressures = 0)
    % secondary-manifold admittance
    Y_rad_xlate = admittance_transfer(TT_sec, Y_rad);
    Y_total_sec = Y_cavity_sec + Nsec*Y_rad_xlate;
    %total
    % Translate passive secondary admittance to primary manifold
    Y_sec_xlate = admittance_transfer(TT_pri_IB*TT_suppr*TT_pri_OB,...
        Y_total_sec);
    Y_total_pri = Y_cavity_pri + (Npri- 1)*Y_sec_xlate;
    %
    % Translate central (primary-total) admittance outbound along
    % primary pipe
    Y_outward_xlate = admittance_transfer(TT_pri_OB*TT_suppr*TT_pri_IB,...
        Y_total_pri);
    %
    % Add (Nsec-1) secondary pipes and secondary cavity
    Y_sec_reference = Y_outward_xlate + Y_cavity_sec +...
        (Nsec-1)*Y_rad_xlate;
    %
    % Compute pressure ratio, p_in/p_sec
    p_rat_sec = [1 0]*TT_rad*TT_sec*[1; Y_sec_reference];
    %
    % Compute pressure ratio, p_sec/p_pri
    p_rat_pri = [1 0]*(TT_pri_OB*TT_suppr*TT_pri_IB)*[1; Y_total_pri];
    %
    % Extra feed pipe to microbarometer: pressure ratio
    p_rat_feed = [1 0]*TT_stub*[1; Y_microbarometer];
    %
    % Compute the pressure ratio, p_pri/p_in
    p_rat(ii) = 1/(p_rat_sec*p_rat_pri*p_rat_feed);
    %
end
% Response component resulting from spatial distribution of inlets
[p_inside, xxbis, yybis, soundvelocity] = response2planewave(freq_vector, ...
    tempC, arrival_theta_radians,...
    Npri, Lpri, Nsec, Lsec, 0);

%
% Combined acoustical and spatial response
p_total_NRS  = Npri * Nsec * p_rat .* p_inside(:, 1);

%=============== Sensor response from MB2005 and MB3
switch firflag
    case'fir'
        TF_sensor_UT                = idc2fap_is(freq_vector,sensor_UT,'fir');
        TF_ref_sensor               = idc2fap_is(freq_vector,ref_sensor,'fir');
        p_totalNRSsensor            = p_total_NRS .* TF_sensor_UT;
        TF_ref_sensor4freqRatio     = idc2fap_is(allfrqsPfilters,ref_sensor,'fir');
    case 'nofir'
        TF_sensor_UT                = idc2fap_is(freq_vector,sensor_UT);
        TF_ref_sensor               = idc2fap_is(freq_vector,ref_sensor);
        p_totalNRSsensor            = p_total_NRS .* TF_sensor_UT;
        TF_ref_sensor4freqRatio     = idc2fap_is(allfrqsPfilters,ref_sensor);
end

