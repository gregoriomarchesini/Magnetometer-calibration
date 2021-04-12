% The script is able to complete a first order Calibration for the 
% 'MIST MOCK MAGNETOMETER TEST'

% PROCEDURE :
% –––––––––
% 
% The magnetometer is exposed to the geomagnetic field which is estimated
% by means of the IGRF model. The three axes of the mgnetometer are aliged
% to to the nord component of the geomagnetic field in sequence. For axis
% the noth componenet of the geomagnetic field is measured by aligning the
% axis with the north direction and taking the measurement both with the
% magnetomet axis 'in' the north direction and and 'opposed' to the north
% direction (south). The procedure gives a total of 6 measurements which
% will be used in order to obtain both 'offset' and 'scale factor' for each
% axis.

% NOTE : remeber all the the measurment are taken with the axes aligned to
% the noth direction so only the north direction of the IGRF model will be
% used for the calibration

% Location of the calibration and IGRF Refernce :
% –––––––––––––––––––––––––––––––––––––––––––––
%
% Location   : Stockholm, KTH University
% Latitude   : 59° 21° 03''
% Longitude  : 18° 04' 23''
% Day        : 18 March 2021  (Thursday)


%% Calibration Model For Each axis 


%    H_sensed_pos  = ( M_real + Offeset)*Scale
%    H_sensed_neg  = (-M_real + Offeset)*Scale

%    H_real   = IGRF megntic field in the direction of evaluation
%    Scale    = Scale factor
%    H_sensed = Measured magnetic field (pos== aligned, neg== opposite)
%    Offset   = Offset in the measurement 

% NOTE : In the following calibration the reference field is the
%        full-intensity geomagnetic field 

%% Main Script

%%  Calibration  Parameters from IGRF Model :

dec_mag = 6.880;                % [deg] Declination ((+)eastward)          
inc_mag = 72.798;               % [deg] Inclination ((+)downward)
M_int   = 51668*10^-9*10^6;     % [microT] 
                           
%% Magnetic field decomposition

M_down  = M_int*sind(inc_mag);                 %(positive downward)
M_est   = M_int*cosd(inc_mag)*sind(dec_mag);
M_north = M_int*cosd(inc_mag)*cosd(dec_mag);

%% MEASURED OUTPUTS : TO BE INSERTED MANUALLY

magx_axis_pos_test =    124.941;    % [microT] microTesla
magx_axis_neg_test =    -101.53 ;   % [microT] microTesla

magy_axis_pos_test =    90.9156;    % [microT] microTesla
magy_axis_neg_test =   -99.2445 ;   % [microT] microTesla

magz_axis_pos_test =    63.3693;    % [microT] microTesla
magz_axis_neg_test =   -155.81;     % [microT] microTesla



%% Determination of the scakle factor 

% NOTE: The scale factor will be considered the same for both directions of
%       the axis. 

Scale_x   = (magx_axis_pos_test-(magx_axis_pos_test+magx_axis_neg_test)/2)/M_int

Scale_y   = (magy_axis_pos_test-(magy_axis_pos_test+magy_axis_neg_test)/2)/M_int

Scale_z   = (magz_axis_pos_test-(magz_axis_pos_test+magz_axis_neg_test)/2)/M_int




%% Offeset Determinbation due to Strong Ferromgnetic Interference

%  NOTE: This offset must be Added (not subtracted) from the 
%         sensed magnetic field 
%  Note: The Offset can be measured indipendently from the scale factor
%        Since no real reference field is required for that



offset_x = (magx_axis_pos_test + magx_axis_neg_test)/2/Scale_x

offset_y = (magy_axis_pos_test + magy_axis_neg_test)/2/Scale_y

offset_z = (magz_axis_pos_test + magz_axis_neg_test)/2/Scale_z





%% Test of the outputs 

M_real_int = (magx_axis_pos_test/Scale_x-offset_x)

save('Calibration_data','offset_x','offset_y','offset_z', ...
     'Scale_x','Scale_y','Scale_z')