%% Returns normalized complex frequency response from IDC response format (for Infrasound)
% on input frequency vector f (Hz)
% and IDC format response file
% by Benoit Doury
%
% Options:
% include_fir = {'fir'|'nofir'} (Default: nofir)
%
% OUT
% TF: the transfer function in its complex form
%
% Example
% >> [TF]=idc2fap_is(f,'I26DE_BDF_RSP_2015134_C','fir');
% >> figure;
% >> semilogx(f, abs(TF));
%

function [TF]=idc2fap_is(f,idc_rsp_file, include_fir)

if ~exist(idc_rsp_file,'file')
    fprintf('this file does not exist\n');
end

%% opening file and initiating variables

fid = fopen (idc_rsp_file, 'r');
tline = fgetl(fid);
jw=1j*f*(2*pi);   %Complex angular frequency vector, in rad/s

%Initiate TF vector
TF = jw./jw;

% Use FIR or not?
if exist('include_fir','var') && strcmp(include_fir,'fir')
    include_fir = 1;
else
    include_fir = 0;
end

stage_id = 0;

% Now going through the file
while ischar(tline)
    if (strfind(tline,'#'))
    else
        if (strfind(tline,'paz'))                     %entering a PAZ stage
            stage_id = stage_id + 1;
            % fprintf('paz, stage %i\n',stage_id);
            
            % 1 - Read the PAZ stage into vectors
            tline = fgetl(fid);
            A0 = textscan(tline,'%f');                %read normalization factor
            tline = fgetl(fid);
            nb_poles = textscan(tline,'%f');          %read number of poles
            p_vector = zeros(nb_poles{1},1);
            for j=1:nb_poles{1},
                tline = fgetl(fid);
                pole = textscan(tline,'%f %f %f %f'); %read poles
                p_vector(j) = pole{1} + 1j*pole{2};
            end
            tline = fgetl(fid);
            nb_zeros = textscan(tline,'%f');
            z_vector = zeros(nb_zeros{1},2);
            for j=1:nb_zeros{1},
                tline = fgetl(fid);
                zero = textscan(tline,'%f %f %f %f'); %read zeroes
                z_vector(j) = zero{1} + 1j*zero{2};
            end
            
            % 2 - Apply PAZ to Transfer Function
            D_TF=jw./jw;                              % Initiate Denominator of TF
            for k=1:nb_poles{1},
                D_TF = ((jw) - p_vector(k)).*D_TF;
            end
            N_TF = jw./jw;                            % Initiate Numerator of TF
            for k=1:nb_zeros{1},
                N_TF = ((jw) - z_vector(k)).*N_TF;
            end
            TF = A0{1}*TF.*N_TF./D_TF;                % Applying PAZ stage to TF
        end
        
        if (strfind(tline,'fir'))                    %entering a FIR stage
            if(include_fir == 1)
                stage_id = stage_id + 1;
                %                 fprintf('fir stage %i\n',stage_id);
                
                %1 - Read the FIR into vector
                tline = fgetl(fid);
                Decimation = textscan(tline,'%f');   %read decimation factor
                tline = fgetl(fid);
                nb_elements = textscan(tline,'%f'); %read number of FIR elements
                fir_vector = zeros(nb_elements{1},1);
                for j=1:nb_elements{1},
                    tline = fgetl(fid);
                    fir_element = textscan(tline,'%f %f');
                    fir_vector(j) = fir_element{1};
                end
                fgetl(fid);                 %getting rid of the trailing 0 in FIR filter.
                
                %2 - Apply FIR to Transfer Function
                delay=grpdelay(fir_vector);
                %                 fprintf('grpdelay %f \n',delay(1)/40);
                fir=freqz(fir_vector,1,f,Decimation{1});
                TF=TF.*fir;
            end
        end
        if (strfind(tline,'fap'))               %entering a FAP stage
            stage_id = stage_id + 1;fprintf('fap stage %i\n',stage_id);
            
            %1 - Read the FAP into vector
            tline = fgetl(fid);
            nb_elements = textscan(tline,'%f'); %reading the number of elements of the FIR filter
            fap_vector = zeros(nb_elements{1},3);
            for j=1:nb_elements{1},
                tline = fgetl(fid);
                fap_element = textscan(tline,'%f %f %f %f %f');
                fap_vector(j,1) = fap_element{1};fap_vector(j,2) = fap_element{2};fap_vector(j,3) = fap_element{3};
            end
            %2 - Apply FAP to Transfer Function
            amp=interp1(fap_vector(:,1),fap_vector(:,2),f);
            phas=interp1(fap_vector(:,1),fap_vector(:,3),f);
            fap=amp.*exp(1j*phas*pi/180);
            TF=TF.*fap;
        end
        %     disp(tline)
    end
    tline = fgetl(fid);
end

fclose (fid);
end



