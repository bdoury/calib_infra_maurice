function [p_total_NRSsensor, p_total_NRS ] = ...
    HCP_acoustical_benoit(freq_vector,...
    sensor_UT, firflag)

%==========================================================================
%              Thermo-viscous acoustical model for a
%             hexagonal close-packed (HCP) pipe system
%
%      This version contains an example system based on the large HCP
%   at the Conrad Observatory (I99).  It is easily modified for other
%   HCP systems.
%
%                  T. Gabrielson (tbg3@psu.edu)
%                  Penn State University
%                  [May 2012]
%
%==========================================================================
% 
% dir ../00benoit/
% addpath ../00benoit/
%================= USER INPUT BLOCK =======================================
%                             (V2_BFH is 18m HCP at I99)
% HCP_aperture = 36;  HCP_order = 3;  %  36m HCP at I99
% HCP_aperture = 18;  HCP_order = 2;  %  18m HCP at I99
%
arrival_elevation_degrees = 30;  %  Elevation angle: horizontal = 0
%    Elevation angle can be a vector.  If it is entered as a vector, then
%  the response is computed and saved for each angle.  [Max 8 angles]
arrival_azimuth_degrees   = 0;  %  Azimuth angle: 0 aligns with
%                                        first primary pipe
%     Can also enter a vector of azimuths for averaging over azimuth:
%  arrival_azimuth_degrees   =  0:2:60;  % To average over the primary
%                                          pipe interval, which is 60 deg
%
%  Pipe specs:
%        numbers and lengths [m] of primary, secondary, and tertiary pipes
% Npri = 6;       Nsec = 3;      Nter = 6;      % number of pipes: 36-m HCP
% L_pri = 10.36;  L_sec = 4.97;  L_ter = 3.08;  % lengths of pipes 36-m HCP
%
% Npri = 6;      Nsec = 2;      Nter = 6;     % number of pipes: 18-m HCP
% L_pri = 3.50;  L_sec = 3.48;  L_ter = 2.07; % lengths of pipes 18-m HCP
%


%====================================================================
Npri  = 4;     Nsec = 12;       Nter = 2;     % number of pipes: 18-m HCP
%====================================================================
%====================================================================
%L_pri = 4.25;  L_sec = 3.00;    L_ter = 0.85;%2.07;  % lengths of pipes 18-m HCP
L_pri = 5.75;  L_sec = 3.05;    L_ter = 0.85/2;%2.07;  % lengths of pipes 18-m HCP
%====================================================================
%====================================================================
% the theoretical model does not work with H1 to H5 in IS26
%     WARNING: change values for 1 to 5
% On other way could be to consider that the "ground truth" 
% can be obtained by averaging of the measureemnts
% but in this case the dip artefact appears.
% real values are:
%=== here L_pri = 5.75;  L_sec = 3.00;    L_ter = 0.85/2;%2.07;  % lengths of pipes 18-m HCP
%
%====================================================================
%====================================================================

%
%  Primary and secondary pipe radii (not diameter) [m]
a_pri = 13e-3/2;%0.824*25.4e-3/2; %12e-3/2;  %  DN12
a_sec = 13e-3/2;%0.824*25.4e-3/2; %12e-3/2;  %0.622*25.4e-3/2;   DN12
a_ter = 13e-3/2;%12e-3/2;  %  DN12
%
%  Length of stub from primary manifold to microbarometer [m]
L_stub = 1.6;%12.24 + 1.16;
a_stub = a_sec;%36e-3/2; % 36-m HCP
% % L_stub = 12.24 + 1.16;  a_stub = a_sec; % 18-m HCP
%
%  Cavity specs (microbarometer, primary, secondary,
%                     and tertiary manifolds)
%
a_microbarometer = 0.0545/2;  %  radius and length to get correct MB3
L_microbarometer = 0.0545;    
% a_ter_cavity = 5/100;   %85e-3/2;  % 85 mm diameter
% L_ter_cavity = 3.8/100; %50e-3;    % 50 mm height
% a_sec_cavity = a_ter_cavity;
% L_sec_cavity = L_ter_cavity;

a_ter_cavity = 13e-3/2;   %85e-3/2;  % 85 mm diameter
L_ter_cavity = 50e-3;     %50e-3;0e-3    % 50 mm height

a_sec_cavity = 53.4e-3;
L_sec_cavity = 112e-3;
%   
if L_stub ~= 0
    L_pri_cavity = 0.02;    % main manifold is a cube of side 4cm
end
%
%   Resonance suppressor specs
%        If either of the 'RS_at_..._manifold' flags are false, the
%          corresponding radius and length inputs are ignored.
RS_at_secondary_manifold = false;
a_rss = 0.975e-3;  L_rss = 20e-3;  %  Radius and length of capillary [m]
%   RS_at_secondary_manifold means that the suppressors are at the
%       outboard ends of the primary pipes.
RS_at_primary_manifold = false;
a_rsp = 0.975e-3;  L_rsp = 20e-3;  %  Radius and length of capillary [m]
%   RS_at_primary_manifold means that a suppressor is at the
%       primary-manifold end of the stub to the microbarometer.
%
%   Temperature [deg C] and ambient pressure [Pa]
tempC = 10;
atmos_press = 88000;  %   88000 Pa is roughly "standard"
%                                  atmospheric pressure at the altitude of
%                                  Conrad Observatory
%
%   Select frequencies for solution
%
%
%     'damping_correction' is an empirical correction factor used in
%  Tmatrix_air.  Set to zero to disable that correction.
damping_correction = 0;  %  normally set to one
%
%
%=============== END USER INPUT BLOCK =====================================
%
%  Required m-functions:
%        Tmatrix_air
%        radiation_admittance
%        Y_cavity
%        Z_capillary_inlet
%        admittance_transfer
%        Inlet_distribution_rosette
%        PW_response
%
%  Miscellaneous pre-calculations
tempK =  tempC + 273;
cc     = 20.05*sqrt(tempK);
MW_air = 28.97;
Rgas   = 8314.3;
rho    = atmos_press*MW_air/(Rgas*tempK);
%
%     The ultimate model result is the ratio of (complex) pressure at the
%   microbarometer to acoustic pressure at the inlet:  p_ratio
p_ratio = zeros(size(freq_vector));  %  Allocate space for pressure ratio
%
%   Set length increment for pipe model based on highest frequency
min_wavelength = cc/freq_vector(end);
max_dx = min_wavelength/1000;
%
%     Set some transfer matrices to the identity matrix.  This permits
%   omitting the feature (by leaving the matrix as the identity matrix or
%   including the feature by changing the values.
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
    Nslice = floor(L_pri/max_dx) + 1;  %  number of "slices"
    TT = Tmatrix_air(freq, a_pri, L_pri/Nslice, tempK,...
        atmos_press, damping_correction); % transfer matrix for pipe slice
    TT_pri = TT^Nslice;  %  transfer matrix for entire length of pipe
    %
    %         Secondary pipe
    Nslice = floor(L_sec/max_dx) + 1;
	TT = Tmatrix_air(freq, a_sec, L_sec/Nslice, tempK,...
        atmos_press,damping_correction);
    TT_sec = TT^Nslice;
    %
    %         Tertiary pipe
    Nslice = floor(L_ter/max_dx) + 1;
    TT = Tmatrix_air(freq, a_ter, L_ter/Nslice, tempK,...
        atmos_press, damping_correction);
    TT_ter = TT^Nslice;
    %
    %         Central-manifold-to-microbarometer pipe (if used)
    if L_stub > 0  %  otherwise, leave matrix as identity matrix
        Nslice = floor(L_stub/max_dx) + 1;
        TT = Tmatrix_air(freq, a_stub, L_stub/Nslice, tempK,...
            atmos_press, damping_correction);
        TT_stub = TT^Nslice;
    end
    %
    %  Radiation admittance and corresponding transfer matrix
    Y_rad  = radiation_admittance(rho, cc, a_ter, freq);
    TT_rad = [1  1/Y_rad;  0  1];
    %
    %  Cavity admittances (assume cylindrical cavities)
    Y_cavity_ter     = Y_cavity(freq, a_ter_cavity, L_ter_cavity,...
        tempK, atmos_press);
    Y_cavity_sec     = Y_cavity(freq, a_sec_cavity, L_sec_cavity,...
        tempK, atmos_press);
    if L_stub == 0
        Y_microbarometer = 0; % microbarometer is primary manifold    
        Y_cavity_pri     = Y_cavity(freq, a_microbarometer, L_microbarometer,...
        tempK, atmos_press);
    else
        Y_microbarometer = Y_cavity(freq, a_microbarometer,...
            L_microbarometer, tempK, atmos_press);    
        Y_cavity_pri     = Y_cavity_cube(freq, L_pri_cavity,...
            tempK, atmos_press);
    end
    %
    %  Resonance suppressors if used
    if RS_at_secondary_manifold  %  otherwise, matrix is identity matrix
        Zrs = Z_capillary_inlet(freq, a_rss, L_rss, tempK, atmos_press);
        TT_RS_sec = [1 Zrs; 0 1];
    end
    if RS_at_primary_manifold  %  otherwise, matrix is identity matrix
        Zrs = Z_capillary_inlet(freq, a_rsp, L_rsp, tempK, atmos_press);
        TT_RS_pri = [1 Zrs; 0 1];
    end
    %
    %     Total "passive" (i.e., all inlet pressures = 0)
    %               tertiary-manifold admittance
    Y_rad_xlate = admittance_transfer(TT_ter, Y_rad);
    Y_total_ter = Y_cavity_ter + Nter*Y_rad_xlate;
    %
    %     Translate passive tertiary admittance to secondary manifold
    Y_ter_xlate = admittance_transfer(TT_sec, Y_total_ter);
    Y_total_sec = Y_cavity_sec + Nsec*Y_ter_xlate;
    %
    %     Translate passive secondary admittance to primary manifold
    Y_sec_xlate = admittance_transfer(TT_pri*TT_RS_sec, Y_total_sec);
    %
    %     Translate microbarometer admittance to primary manifold
    Y_micro_xlate = admittance_transfer(TT_RS_pri*TT_stub,...
        Y_microbarometer);
    Y_total_pri = Y_cavity_pri + Y_micro_xlate + (Npri - 1)*Y_sec_xlate;
    %
    %     Translate central (primary-total) admittance outbound along
    %               primary pipe
    Y_outpri_xlate = admittance_transfer(TT_RS_sec*TT_pri, Y_total_pri);
    %
    %     Add Nsec - 1 secondary pipes and secondary cavity
    Y_sec_reference = Y_outpri_xlate + Y_cavity_sec +...
        (Nsec - 1)*Y_ter_xlate;
    %
    %     Translate Y_sec_reference outbound along secondary pipe
    Y_outsec_xlate = admittance_transfer(TT_sec, Y_sec_reference);
    %
    %     Add Nter - 1 tertiary pipes and tertiary cavity
    Y_ter_reference = Y_outsec_xlate + Y_cavity_ter +...
        (Nter - 1)*Y_rad_xlate;
    %
    %================================================
    %
    %  Compute pressure ratio, p_inlet/p_ter
    p_ratio_ter  = [1 0]*TT_rad*TT_ter*[1; Y_ter_reference];
    %
    %  Compute pressure ratio, p_ter/p_sec
    p_ratio_sec  = [1 0]*TT_sec*[1; Y_sec_reference];
    %
    %  Compute pressure ratio, p_sec/p_pri
    p_ratio_pri  = [1 0]*TT_RS_sec*TT_pri*[1; Y_total_pri];
    %
    %  Extra feed pipe to microbarometer: pressure ratio
    p_ratio_feed = [1 0]*TT_RS_pri*TT_stub*[1; Y_microbarometer];
    %
    %  Compute the pressure ratio, p_pri/p_in
    p_ratio(ii) = 1/(p_ratio_ter*p_ratio_sec*p_ratio_pri*p_ratio_feed);
    %
end
%
%   Sum over all inlets (equivalent to vertical incidence)
p_ratio = p_ratio*Npri*Nsec*Nter;
%
%   End of primary computational loop (one pass per frequency)
%
%      Overall response by combining the results of all inlet excitations
%   (superposition)
%
%        elevation angle (zero is horizontal)
arrival_elevation_radians = arrival_elevation_degrees*pi/180;
%        azimuth angle (zero aligns with first primary pipe)
arrival_azimuth_radians   = arrival_azimuth_degrees*pi/180;
%
%   Generation of overall response by superposition
%
%       Calculate X,Y coordinates for all inlets
% [xx, yy] = Inlet_distribution_HCP(HCP_aperture, HCP_order);

[xx, yy] = Inlet_distribution_T_rosette(Npri, L_pri,...
    Nsec, L_sec);

%
%       Find plane-wave response based on the distribution of inlets.
%          Note: if more than one elevation angle is specified, the
%        response is calculated and stored for each elevation angle; if
%        more than one azimuth angle is specified, the response is
%        calculated as the average response over those azimuth angles.
%
p_inside = PW_response(freq_vector, cc,...
    arrival_elevation_radians, arrival_azimuth_radians, xx, yy);
%   Combined acosutical and spatial response
p_total_NRS  = p_ratio .* p_inside;

%=============== Sensor response from MB2005 and MB3
switch firflag
    case'fir'
        TF_sensor_UT                = idc2fap_is(freq_vector,sensor_UT,'fir');
        p_total_NRSsensor            = p_total_NRS .* TF_sensor_UT;
    case 'nofir'
        TF_sensor_UT                = idc2fap_is(freq_vector,sensor_UT);
        p_total_NRSsensor            = p_total_NRS .* TF_sensor_UT;
end
