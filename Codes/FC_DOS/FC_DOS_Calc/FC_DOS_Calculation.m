%Code to calculate the lineshape function and Franck-Condon-weighted density of states (FC-DOS, F(wab))
%Based on May and Kuhn 2000, Eqs 5.1.39 through 5.1.43
    % Assumptions: 1) normal modes of ground state are used for all excited states
    %              2) harmonic approximation
    %INPUT - 1) transition energy levels (energy_pop), 2)Huang-Rhys factors (output_HR_matlab), 3) phonon modes (modes_cm)
    %--Can change output based on smearing parameter for temperature (temp_smearing)
%OUTPUT - 1) DeltaG(t), 2) DeltaG(0), 3) Lineshape function 4) F(wab) matrix
%OUTPUT - Key output is fcf_weighted_dos_values.txt, which is used as input for relaxation dynamics code

%NOTES: 
%1) You must confirm that the sampling frequency and number of time steps
%used for the inverse Fourier transform give positive values for the real
%part of the transform ('ft_exp_full_Gt').  In general, a sampling
%frequency of 10 Hz and total number of time steps of 10000 to 100000
%should be sufficient
%2) The phonon modes used in the Fourier Transform must have only a small number of
%significant figures - otherwise a huge sampling frequency and number of
%time steps is required for the Fourier transform to provide positive
%values as discussed in 1).  Therefore, the code converts the input modes
%to units of 1/fs and then rounds to a certain number of sig figs (3 by default)
%3) Adjust smearing parameters for different finite temperatures

%--------Constants-----------%
hbar = 0.65821; % Planck's constant in eV*fs
hnobar = 4.135; % Planck's constant without bar in eV*fs
%-----------INPUT-------------%

%Energy levels
e_min = 1; % First excited state index
e_max = 7; % Last excited state index
energy_pop = load('energy_pop')*[0 1 0]'; % System energy levels in eV

%Huang-Rhys factors for each phonon mode (row) and transition energy (column)
HR_factors = load('output_HR_matlab'); % dimensionless

%Normal modes from ground state calculation in same order as output in Gaussian
modes_cm = load('modes_cm'); % modes in cm^-1 from phonon mode calculation
nmodes = length(modes_cm); % total number of modes
modes_ev = zeros(nmodes,1); % modes in eV
modes_fs = zeros(nmodes,1); %modes in 2*pi/fs
for i=1:nmodes;
    modes_ev(i) = modes_cm(i)/8065.73; % modes in eV, converted from gaussian cm-1
    modes_fs(i) = modes_ev(i)/hbar; %modes in 2*pi/fs
    modes_nu(i) = modes_fs(i)/(2*pi); %modes in 1/fs
end;

%The modes must be read in with a limited number of significant figures,
%otherwise an enormous number of timesteps will be required to capture
%the exact period of each frequency.  Therefore, he we round modes_nu to
%a certain number of sig figs to be used as input.  Not doing this will
%lead to negative values in the fourier transform.  Always
%check your data to make sure the real part (real(ifft_exp_full_Gt)) is
%positive!  In general, 3 sig figs works for about 10000 timesteps.
for i=1:nmodes;
    modes_nu(i) = round(modes_nu(i),3,'decimal');
end;

% Time step for Fourier transform:
num_ts = 10000; % converge this value and make sure it gives positive values of real part of Fourier transform!
f_s = 10; % sampling frequency in Hz - this will determine energy intervals along which F(wab) is measured

%Set up sampling frequency and length of FT output
ts_int = 1/f_s;  % Sampling frequency
ft_max_energy = 5; % Limit for cutting off Fourier transform result for Gaussian smearing and fittin g(in eV)

% Create time vector for Fourier transform
time_steps = zeros(num_ts,1);
for i=1:num_ts;
    time_steps(i) = ts_int*i;
end;

%Smearing
temp_smearing = 0.0257; % temperature smearing for B-E distribution in eV

%Optional Lorentzian smearing based on dephasing time - choosing 0 for both will provide unsmeared F(w):
smearing_l = 0; % 1 for yes; 0 for no
l_smearing = 0.0257; % Lorentzian smearing in eV
%Optional Gaussian smearing for thermal effects:
smearing_g = 0; % 1 for yes; 0 for no
g_smearing = 0.0257; %Gaussian smearing in eV

%-----------END INPUT-------------% 

%Calculate matrix of transition energies for all band combinations in eV!
%NOTE: take absolute value to keep all transition energies positive - this
%is because we only calculate downward transitions - upward ones are
%calculated using detailed balance in relaxation dynamics code
num_elevels = e_max - e_min + 1;
band_comb = (num_elevels)*(num_elevels - 1)/2;
band_diff = zeros(num_elevels,num_elevels);
band_diff_one_row = zeros(1,band_comb); %same value but in one column
band=1;
for i=1:num_elevels;
    for l=(i+1):num_elevels;
        band_diff(i,l) = (abs(energy_pop(l) - energy_pop(i))); %eV
        band_diff(l,i) = band_diff(i,l);
        band_diff_one_row(1,band) = band_diff(i,l);
        band = band+1;
    end;
end;

%Calculate Bose-Einstein distribution for all phonon modes
BE_dist = zeros(nmodes,1);
for i=1:nmodes;
    BE_dist(i,1) = 1/(exp(modes_ev(i)/temp_smearing)-1);
end;

%Create integrand for Fourier transform (see May and Kuhn)
full_G = zeros(num_ts,band_comb);
exp_full_G = zeros(num_ts,band_comb);
for i=1:band_comb;
        for t=1:num_ts;
            full_G(t,i) = 0.0;
            for l=1:nmodes;
                full_G(t,i) = full_G(t,i) + HR_factors(l,i)*((exp(-1j*(2*pi*modes_nu(l))*time_steps(t))-1)*(1 + BE_dist(l)) + (exp(1j*(2*pi*modes_nu(l))*time_steps(t))-1)*BE_dist(l));
            end;
            exp_full_G(t,i) = exp(full_G(t,i));
        end;
end;

%Set up x-axes of different units
x_axis_fs = zeros(num_ts,1); % x axis in 1/fs
x_axis_ev = zeros(num_ts,1); % x-axis in eV
x_axis_cm = zeros(num_ts,1); % x-axis in cm^-1
x_axis_fs(1) = 0.0;
for t=2:num_ts;
    x_axis_fs(t) = x_axis_fs(t-1) + f_s/num_ts; %x-axis in 1/fs
end;
for t=1:num_ts;
    x_axis_ev(t) = x_axis_fs(t)*1000000000000000*2*3.14159*0.0000000000000006582;
    x_axis_cm(t) = x_axis_fs(t)*1000000000000000*2*3.14159*0.0000000000000006582/0.0001239;
end;

%Compute Fourier transform: multiply by timestep for correct units and scaling
%Scale by the number of timesteps and the time-step interval to get correct amplitude (matlab divides ifft by num_ts implicitly)
%This scaling ensures that Parseval's Theorem will hold, which is probably the best systematic way to scale Fourier transforms!
ifft_exp_full_Gt = ifft(exp_full_G)*ts_int*num_ts;
ft_exp_full_Gt = real(ifft_exp_full_Gt);
ft_back_exp_full_Gt = fft(ifft_exp_full_Gt)*f_s/num_ts;

ft_expDeltaGt = ft_exp_full_Gt;

%Multiply Fourier transform by 1/hbar prefactor to get units of 1/eV:
ft_expDeltaGt_evfs = zeros(num_ts,band_comb);
for i=1:band_comb;
        ft_expDeltaGt_evfs(:,i) = ft_expDeltaGt(:,i)/hbar;
end;

%Create time series only up to cutoff used for Fourier Transform
%(Fourier transform provides energy range significantly greater than typical transition energies):
num_cutoff = round(ft_max_energy/(f_s/num_ts*1000000000000000*2*3.14159*0.0000000000000006582)); % max timestep for cutoff
x_axis_fs_cut = zeros(num_cutoff,1); % x axis in 1/fs up to cutoff value
x_axis_ev_cut = zeros(num_cutoff,1); % x-axis in eV up to cutoff value
x_axis_cm_cut = zeros(num_cutoff,1); % x-axis in cm-1 up to cutoff value

x_axis_fs_cut(1) = 0.0;
for t=2:num_cutoff;
    x_axis_fs_cut(t) = x_axis_fs_cut(t-1) + f_s/num_ts; %x-axis in 1/fs
end;
for t=1:num_cutoff;
    x_axis_ev_cut(t) = x_axis_fs_cut(t)*1000000000000000*2*3.14159*0.0000000000000006582;
    x_axis_cm_cut(t) = x_axis_fs_cut(t)*1000000000000000*2*3.14159*0.0000000000000006582/0.0001239;
end;

ft_expDeltaGt_evfs_cutoff = zeros(num_cutoff,band_comb);
for t = 1:num_cutoff;
    ft_expDeltaGt_evfs_cutoff(t,:) = ft_expDeltaGt_evfs(t,:);
end;

ft_expDeltaGt_evfs_cut_FINAL = zeros(num_cutoff,band_comb);
ft_expDeltaGt_evfs_cut_FINAL = ft_expDeltaGt_evfs_cutoff;

%---OPTIONAL SMEARING (optionally set with 'smearing_l' and 'smearing_g' in input section)-----%
ft_expDeltaGt_evfs_cut_gauss = zeros(num_cutoff,band_comb);
ft_expDeltaGt_evfs_cut_gauss_norm = zeros(num_cutoff,band_comb);
ft_expDeltaGt_evfs_cut_lorentz = zeros(num_cutoff,band_comb);
ft_expDeltaGt_evfs_cut_lorentz_norm = zeros(num_cutoff,band_comb);

%A) Lorentzian smearing ('l_smearing' input parameter determines broadening)
if (smearing_l==1)
    lorentz_sum = zeros(num_cutoff,1);
    for i=1:band_comb;
        for j=1:num_cutoff;
            lorentz_sum(j,1) = 0.0;
            for k=1:num_cutoff;
                    ft_expDeltaGt_evfs_cut_lorentz(j,i) = ft_expDeltaGt_evfs_cut_lorentz(j,i) + ft_expDeltaGt_evfs_cutoff(k,i)*((l_smearing/2)/((l_smearing/2)^2+((x_axis_ev_cut(j)-x_axis_ev_cut(k))^2)));
                    lorentz_sum(j,1) = lorentz_sum(j,1) + ((l_smearing/2)/((l_smearing/2)^2+((x_axis_ev_cut(j)-x_axis_ev_cut(k))^2)));
            end
            ft_expDeltaGt_evfs_cut_lorentz_norm(j,i) = ft_expDeltaGt_evfs_cut_lorentz(j,i)/lorentz_sum(j,1);
        end
    end
    ft_expDeltaGt_evfs_cut_FINAL = ft_expDeltaGt_evfs_cut_lorentz_norm;
end
            
%B) Gaussian smearing ('g_smearing' input parameter determines broadening)
if (smearing_g==1)
    gaussian_sum = zeros(num_cutoff,1);
    for i=1:band_comb;
        for j=1:num_cutoff;
            gaussian_sum(j,1) = 0.0;
            for k=1:num_cutoff;
                    ft_expDeltaGt_evfs_cut_gauss(j,i) = ft_expDeltaGt_evfs_cut_gauss(j,i) + ft_expDeltaGt_evfs_cutoff(k,i)*exp(-((x_axis_ev_cut(j)-x_axis_ev_cut(k))^2)/(2*g_smearing*g_smearing));
                    gaussian_sum(j,1) = gaussian_sum(j,1) + exp(-((x_axis_ev_cut(j)-x_axis_ev_cut(k))^2)/(2*g_smearing*g_smearing));
            end
            ft_expDeltaGt_evfs_cut_gauss_norm(j,i) = ft_expDeltaGt_evfs_cut_gauss(j,i)/gaussian_sum(j,1);
        end
    end
    ft_expDeltaGt_evfs_cut_FINAL = ft_expDeltaGt_evfs_cut_gauss_norm;
end

%----END OPTIONAL SMEARING SECTION----%

%Use cubic spline interpolation to calculate value of lineshape at
%transition energy for each band combination
ls_spline_value_one_row = zeros(1,band_comb); %Variable providing lineshape function's value at transition energy
for i=1:band_comb;
    ls_spline_value_one_row(1,i) = interp1(x_axis_ev_cut(:,1),ft_expDeltaGt_evfs_cut_FINAL(:,i),band_diff_one_row(i),'spline');
end;

%Convert single row of lineshape values at transition energies
%(band_diff_one_row) to matrix to input into relaxation dynamics
ls_spline_value_mat = zeros(num_elevels,num_elevels);
band = 1;
for i=1:num_elevels;
    for l=(i+1):num_elevels;
        ls_spline_value_mat(i,l) = ls_spline_value_one_row(1,band);
        ls_spline_value_mat(l,i) = ls_spline_value_mat(i,l);
        band = band+1;
    end;
end;

%Write important output to file:

%x_axis_ev
dlmwrite('x_axis_ev.txt',x_axis_ev,'precision',5,'delimiter',' ')

%x_axis_cm
dlmwrite('x_axis_cm.txt',x_axis_cm,'precision',5,'delimiter',' ')

%ft_expDeltaGt_evfs: this if F(w) - value at F(wab) is matrix element for that transition entering into relaxation dynamics
dlmwrite('lineshapes.txt',ft_expDeltaGt_evfs_cut_FINAL,'precision',5,'delimiter',' ')

%ls_spline_value_mat - this output matrix is F(wab) - used as input for relaxation dynamics
dlmwrite('fcf_weighted_dos_values.txt',ls_spline_value_mat,'precision',5,'delimiter',' ')
