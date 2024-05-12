%% LAB 4: RADON TRANSFORM

F = phantom();
%imshow(F);
%disp(size(F));

example = 5; 

%%% the parameter example determines the 
%%% functionality of the code:
%   For example = 
%   1:  MATLAB iradon function is applied for different angle ranges
%   2:  MATLAB iradon function is applied for different filters
%   3:  CUSTOM inverse_radon applied function for different filters and
%       ranges 
%   4:  Alternative CUSTOM inverse_radon function with the linear custom
%       interpolation


if example == 1
    % MATLAB iradon function for different angle ranges

    theta1 = 0:10:170; [R1,xp1] = radon(F,theta1);
    theta2 = 0:5:175;  [R2,xp2] = radon(F,theta2);
    theta3 = 0:2:178;  [R3,xp3] = radon(F,theta3);
    
    figure
    imagesc(theta3,xp3,R3) 
    colormap(hot); colorbar
    xlabel('\theta'); ylabel('x\prime')
    title("Radon Transform of Head Phantom Using 90 Projections")
    
    I1 = iradon(R1,10);
    I2 = iradon(R2,5);
    I3 = iradon(R3,2);
    
    montage({I1,I2,I3},Size=[1 3])

elseif example == 2
    % MATLAB iradon function for different filters

    step_angle = 5;
    theta = 0:step_angle:179;
    [R0, x0] = radon(F, theta);
    figure
    imagesc(theta,x0,R0) 
    colormap(hot); colorbar
    xlabel('\theta'); ylabel('x\prime')
    title("Radon Transform of Head Phantom Using 90 Projections")
    
    I1 = iradon(R0,step_angle, "linear", "Ram-Lak");
    I2 = iradon(R0,step_angle, "linear", "Shepp-Logan");
    I3 = iradon(R0,step_angle, "linear", "Cosine");
    
    montage({I1,I2,I3},Size=[1 3])


elseif example == 3
    % CUSTOM inverse_radon function 

    % Parameters:

    pi_angle = 179;
    p = 367; % npos
    q = 183;
    
    step_angle = 0.1; % dtheta
    step_grid = 1; % step
    
    thetas = 0:step_angle:pi_angle;

    ph = phantom();
    Rph = radon(ph,thetas);
        
    inverse_Rph1 = inverse_radon(Rph, q, thetas, step_grid, 'RamLak_KS');
    inverse_Rph2 = inverse_radon(Rph, q, thetas, step_grid, 'Ram-Lak');
    inverse_Rph3 = inverse_radon(Rph, q, thetas, step_grid, 'RamLak_KS');
    %inverse_Rph3 = inverse_radon(Rph, q, step_angle, step_grid, 'SheppLogan');
    
    montage({inverse_Rph1,inverse_Rph2, inverse_Rph3},Size=[1 3])

elseif example == 4
    % Alternative custom inverse_radon function without the interpolation
    % Parameters:

    ph = phantom();
    Rph = radon(ph,0:179);
        
    inverse_Rph1 = iradon_custom_(Rph,0:179,'Ram-Lak');
    inverse_Rph2 = iradon_custom_(Rph,0:179,'Shepp-Logan');
    inverse_Rph3 = iradon_custom_(Rph,0:179,'Cosine');
    
    montage({inverse_Rph1,inverse_Rph2, inverse_Rph3},Size=[1 3])

elseif example == 5
    % Alternative custom inverse_radon function without the interpolation
    % Parameters:

    ph = phantom();
    Rph = radon(ph,0:179);
        
    inverse_Rph1 = iradon_custom_(Rph,0:179,'Ram-Lak');
    inverse_Rph2 = iradon_custom_(Rph,0:178,'Shepp-Logan');
    inverse_Rph3 = iradon_custom_(Rph,0:2:178,'Cosine');
    
    montage({inverse_Rph1,inverse_Rph2, inverse_Rph3},Size=[1 3])

end


%% FUNCTION DEFINITION 

function [inverse] = inverse_radon(Rph, q, thetas, step_grid, filter)
%   Computes inverse Radon transform:
%   reconstructs the image I from projection 
%   data in the 2-D array R.  
%   
%   thetas describes the angles (in degrees) at which the projections    
%   were taken.  
%   
%   FILTER specifies the filter to use for frequency domain filtering.  
%   FILTER is a string that specifies any of the following standard 
%   filters:
% 
%   'Ram-Lak'     The cropped Ram-Lak or ramp filter (default).  The    
%                 frequency response of this filter is |f|.  Because 
%                 this filter is sensitive to noise in the projections, 
%                 one of the filters listed below may be preferable.   
%   'Shepp-Logan' The Shepp-Logan filter multiplies the Ram-Lak
%                 filter by a sinc function.
%   RamLak_KS'     Modification on the definition of the Ram-Lak filter as
%                  stated by Kak and Shakey, see the definition of the
%                  functions bellow

    p = size(Rph,1);

    rho2 = q; % max rho
    rho1 = -q; % min rho;

    angles = deg2rad(thetas);
    nproj = length(angles);
    step_angle = 0.1;
    %step_angle = nproj/179;
    T = 1;  % sample spacing
    
 
    % Definition of the filter
    h = 0;

    if (strcmp(filter,'Ram-Lak')) 
        h = ramLakFilter(p, T);        
    end

    if (strcmp(filter, 'RamLak_KS'))
        h = ramLak_KS(p, T);
    end 

    if (strcmp(filter, 'Shepp-Logan'))
        h = sheppLoganFilter(p, T);
    end 

    % STEP 1: discrete convolution
   
    for iproj = 1:nproj
        Rph(:, iproj) = conv(Rph(:,iproj), h, 'same');
    end
    filteredProjections = Rph;


    % grid for reconstructed image
    x1 = rho1:step_grid:rho2;
    y1 = rho1:step_grid:rho2;
    [y, x] = ndgrid(-y1, x1);

    % positions made to correspond to N point fft
    drho = (rho2 - rho1) / (p-1);
    rho  = rho1 + (0:(p-1)) * drho;
    positions = zeros(nproj, length(rho));
    for i = 1:nproj
        positions(i, :) = rho;
    end

    % STEP 2: discrete backprojection using a linear interpolation
    % display the image through backprojection
    fdata = zeros(size(x));
    for iproj = 1:nproj
        theta = angles(iproj);
        rho1 = x*cos(theta) + y*sin(theta); % rotate coordinate system by theta
        
        r = positions(iproj,:);
    
        fdata1 = filteredProjections(1:p,iproj); % filtered projections
        
        fdata2 = interp1(r, fdata1, rho1, 'linear', 0); 
        fdata = fdata + deg2rad(step_angle) * fdata2; 
    end

    inverse = fdata;

end 


function [img,H] = iradon_custom_(Rph, thetas, filter)
%   Computes inverse Radon transform:
%   reconstructs the image I from projection 
%   data in the 2-D array R.  
%   
%   thetas describes the angles (in degrees) at which the projections    
%   were taken.  
%   
%   FILTER specifies the filter to use for frequency domain filtering.  
%   FILTER is a string that specifies any of the following standard 
%   filters:
% 
%   'Ram-Lak'     The cropped Ram-Lak or ramp filter (default).  The    
%                 frequency response of this filter is |f|.  Because 
%                 this filter is sensitive to noise in the projections, 
%                 one of the filters listed below may be preferable.   
%   'Shepp-Logan' The Shepp-Logan filter multiplies the Ram-Lak
%                 filter by a sinc function.
%   'Cosine'      The cosine filter multiplies the Ram-Lak filter 
%                 by a cosine function.  

%   References: 
%      A. C. Kak, Malcolm Slaney, "Principles of Computerized Tomographic
%      Imaging", IEEE Press 1988.

    p = size(Rph, 1);
    N = 2*floor(size(Rph,1)/(2*sqrt(2)));
    theta = deg2rad(thetas);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition of the Filter used:

    order = max(64,2^nextpow2(2*p));
    filt = 2*( 0:(order/2) )./order;
    w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

    if (strcmp(filter,'Ram-Lak')) 
        % Do nothing     
    elseif (strcmp(filter, 'Shepp-Logan'))
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2))./(w(2:end)/(2)));
    elseif (strcmp(filter, 'Cosine'))
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2)); 
    end 
    filt(w>pi) = 0;                      % Crop the frequency response
    filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter

    H = filt;
    Rph(length(H),1)=0;  % Zero pad projections 
    Rph = fft(Rph);    % R holds fft of projections


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: discrete convolution
    
    for i = 1:size(Rph,2)
       Rph(:,i) = Rph(:,i).*H; % frequency domain filtering
    end
    
    Rph = real(ifft(Rph));     % filtered projections
    Rph(p+1:end,:) = [];   % Truncate the filtered projections
    img = zeros(N);        % Allocate memory for the image.
    
    % Define the x & y axes for the reconstructed image so that the origin
    % (center) is in the spot which RADON would choose.
    xax = (1:N)-ceil(N/2);
    x = repmat(xax, N, 1);    % x coordinates, the y coordinates are rot90(x)
    y = rot90(x);
    
    costheta = cos(theta);
    sintheta = sin(theta);
    ctrIdx = ceil(p/2);     % index of the center of the projections

    % Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
    % quantity is greater than the length of the projections
    imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.
    if size(Rph,1) < imgDiag 
       rz = imgDiag - size(Rph,1);  % how many rows of zeros
       Rph = [zeros(ceil(rz/2),size(Rph,2)); Rph; zeros(floor(rz/2),size(Rph,2))];
       ctrIdx = ctrIdx+ceil(rz/2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: discrete backprojection using linear interpolation
    for i=1:length(theta)  
      proj = Rph(:,i);
      t = x.*costheta(i) + y.*sintheta(i); 
      a = floor(t);  
      img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
    end
    
    img = img*pi/(2*length(theta));
end



function [inverse] = inverse_radon_2(Rph, q, step_angle, step_grid, filter)
% Modification of the first custom radon function that does not provide
% expected results
    p = size(Rph, 1);

    rho2 = q; % max rho
    rho1 = -q; % min rho;

    angles = deg2rad(0:step_angle:179);
    nproj = length(angles);

    T = 1;  % sample spacing
    
    % Definition of the filter
    h = 0;

    if (strcmp(filter,'RamLak')) 
        h = ramLakFilter(p, T);
         
    end

    if (strcmp(filter, 'RamLak_KS'))
        h = ramLak_KS(p, T);
    end 

    if (strcmp(filter, 'SheppLogan'))
        h = sheppLoganFilter(p, T);
    end 

    % STEP 1: discrete convolution
   
    for iproj = 1:nproj
        Rph(:, iproj) = conv(Rph(:,iproj), h, 'same');
    end
    filteredProjections = Rph;


    % grid for reconstructed image
    x1 = rho1:step_grid:rho2;
    y1 = rho1:step_grid:rho2;
    [y, x] = ndgrid(-y1, x1);

    % positions made to correspond to N point fft
    drho = (rho2 - rho1) / (p-1);
    rho  = rho1 + (0:(p-1)) * drho;
    positions = zeros(nproj, length(rho));
    for i = 1:nproj
        positions(i, :) = rho;
    end

    % STEP 2: discrete backprojection using linear interpolation
    inverse = zeros(size(x));
    for iproj = 1:nproj
        theta = angles(iproj);
        rho1 = x*cos(theta) + y*sin(theta); % rotate coordinate system by theta
    
        for ix = 1:size(x, 1)
            for iy = 1:size(x, 2)
                %x_val = x(ix, iy);
                x_val = rho1(ix, iy);
                k = floor(x_val * theta / drho);
                %disp(k);
                eta = (x_val * theta / drho) - k;
    
                if k < p && k > 0
                    fdata1 = filteredProjections(:, iproj); % filtered projections
                    fdata = (1 - eta) * fdata1(k) + eta * fdata1(k+1);
                    inverse(ix, iy) = inverse(ix, iy) + (2 * pi / p) * fdata;
                end
            end
        end
    end
end


%% RamLak filter function

function h = ramLakFilter(p, T)
    f = -floor(p/2):T:floor(p/2); % Frequency axis
    omega = 2 * pi * f;
    u = sinc(omega) - 0.5 * (sinc(omega/2)).^2;
    phi_hat = zeros(size(f));
    phi_hat(abs(f) <= p/2) = 1;
    v_Omega = phi_hat .* (omega.^2) .* u;
    h = v_Omega;
    %h = ifftshift(v_Omega);  % Shift the filter coefficients for convolution
end


%% RamLak filter KS
function h = ramLak_KS(p, T)
% build filter h according to Kak and Shakey 
    h = -floor(p/2):T:floor(p/2);
    
    for i = 1:length(h)
        if (mod(h(i),2) == 1) 
            h(i) = -1./(h(i)^2 * pi^2 * T^2);
        else
            h(i) = 0;
        end    
        h(ceil(p/2)) = 1/(4*T^2);
    end
end

%% Shepp-Logan Filter

function h = sheppLoganFilter(p, T)
    f = -floor(p/2):T:floor(p/2);  % Frequency axis
    omega = 2 * pi * f;
    u = zeros(size(omega));
    
    % Compute u(s) for non-zero values
    nonzero_idx = abs(omega) ~= pi/2;
    u(nonzero_idx) = (pi/2 - omega(nonzero_idx).*sin(omega(nonzero_idx))) ./ ...
        ((pi/2)^2 - omega(nonzero_idx).^2);
    
    % Compute u(s) for zero values (avoid division by zero)
    zero_idx = abs(omega) == pi/2;
    u(zero_idx) = 1/pi;
    
    phi_hat = sinc(f * pi/2);
    v_Omega = (2 * omega.^2 / pi) .* u .* phi_hat;
    %h = v_Omega;
    h = ifftshift(v_Omega);  % Shift the filter coefficients for convolution
end



