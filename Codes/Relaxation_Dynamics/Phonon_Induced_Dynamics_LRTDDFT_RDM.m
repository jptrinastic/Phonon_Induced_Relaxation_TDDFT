%Code to calculate the phonon-induced relaxation dynamics from an initially excited electron-hole pair within LR-TDDFT
%Input files: 1) energy_pop: contains three columns of state index, energy, and occupancy, respectively
%             2) finalmat: nonadiabatic electronic coupling matrix (V_alpha,beta), calculated using Electronic_NAC codes
%             3) fcf_weighted_dos_values: Franck-Condon-weighted density of states matrix, calculated using FCDOS codes
%Output: 1) PPe: population of each excited state as a function of time

%INPUT--------------------------------------------------------
%-----------MAIN INPUT PARAMETERS--------------:
pi=3.14159; %pi constant

zmax=300; %maximum z distance used for spatial relaxation dynamics (not used in current version)
Z=1:zmax;Zz=Z/zmax*10.12; %units of Z for spatial relaxation visualization (not used in current versioN)

Temperature=300;% temperature in K
kT=0.025/300;%eV = 300 K

energy_pop=load('energy_pop')*[0 1 0]'; %Input TDDFT ENERGY LEVELS
forME=load('finalmat'); %Input SQUARED COUPLING MATRIX from C++ code using TDDFT input
Vab2=forME;

field=load('fcf_weighted_dos_values.txt'); %Input FCF-weighted DOS

forME=2*pi*(abs(forME)).*field./0.658218*1000; % RELAXATION RATE MATRIX in 1/ps

HOMO=0; %Highest Occupied Molecular Orbital taken from energy_pop
LUMO=HOMO+1; %Lowest Unoccupied Molecular Orbita taken from energy_pop

Omin=1; %Minimum band
Omax=7; %Maximum band
inib=5; %Initial excitation, measured from 1st excited state

%-----------END MAIN INPUT PARAMETERS------------

%-----------OUTPUT PARAMETERS----------------------------------

% Define number of timesteps
timeSTEPS=600; %timesteps for dynamics
tgrid=1:timeSTEPS;
%Choose one for time-step spacing
    %time=tgrid/200*3;
    %time=(exp(-7+12.298317*tgrid/(max(tgrid))));     % 0.92fs ~ 200 ps
    %time=(exp(-7+12.70378*tgrid/(max(tgrid))));     % 0.92fs ~ 300 ps
    %time=(exp(-7+13.907755*tgrid/(max(tgrid))));     % 0.92fs ~ 1000 ps
    time=(exp(-7+16.2103404*tgrid/(max(tgrid))));    % 0.92fs ~ 10000 ps
    %time=(exp(-7+18.5131014*tgrid/(max(tgrid))));    % 0.92fs ~ 100000 ps

Egrid=-0.1:.005:3.5; %y-axis for energy visualization
de=.005;
slice=zeros(size(Egrid));
Width=.0257; %smearing for energy visualization

%---------END OUTPUT PARAMETERS-------------------------------

%END INPUT----------------------------------------------------

%Create state vector
Ee=energy_pop;     %conduction band

%Create rate matrix (Re)
Re=forME;  %V_ij
Re_pre = Re; %Initial rate matrices for reference

%Create Redfield population transfer matrices :
Ge=zeros(size(Ee,1),size(Ee,1));

%--------ELECTRON-HOLE REDFIELD POPULATION TRANSFER MATRIX INITIALIZATION:
%Enforce detailed balance by multipling higher energy transitions by exponential factor:
for i=1:Omax-Omin+1;
for j=1:Omax-Omin+1;
    if Ee(j)>Ee(i);
        Re(i,j)=Re(i,j)*exp(-(Ee(j)-Ee(i))/0.0257);
    end;
end;
end;
    Re_db = Re; % Matrix to see effect of detailed balance

%Fill Redfield population transfer matrices based on (3.8.23) in May and Kuhn
for i=1:Omax-Omin+1;
    for j=1:Omax-Omin+1;
        Ge(i,i)=Ge(i,i)-Re(i,j);
    end;
end;
for i=1:Omax-Omin+1;
    for j=i+1:Omax-Omin+1;
        Ge(i,j) = Re(j,i);
        Ge(j,i) = Re(i,j);
    end;
end;

%Create final electron Redfield population transfer matrix for diagonalization:
Re_final = Ge;

%-----DYNAMICS---------------------------------------------------

%Set up initial conditions for electron-hole dynamics
PPe=zeros(timeSTEPS,Omax-Omin+1)';PPe(inib,1)=1;

tmp=PPe;
s=zeros(size((Egrid)'));

%Preparation for standard diagonalization:
psi0_e=PPe(:,1);
Kmax_e=Omax-Omin+1;
t_e=time;

%Standard diagonalization
[U_e,D_e]=eig(Re_final); % Diagonalize Re matrix to get D_e eigenvalues and U_e eigenvector matrix
c_e=U_e\psi0_e;     % inv(U_e)*psi0_e as transformed initial excited state
c_e1=c_e;
FREQ_e=real(diag(D_e));    % Place eigenvalues in vector of size (:,1)
psi_e=c_e(1)*U_e*[1  zeros(1,Kmax_e-1)]'*exp(D_e(1,1)*t_e); % Dynamics
psi_e_mid = psi_e;
for k=2:Kmax_e;
    vec_e=[zeros(1,k-1)  1  zeros(1,Kmax_e-k)]';
psi_e=psi_e+c_e(k)*U_e*vec_e*exp(FREQ_e(k)*t_e);    
end;

psi_e_mid_2 = psi_e;

%Enforce self-consistent normalization
for k=1:10;
    for i=1:Kmax_e;
        psi_e(i,:)=real(psi_e(i,:).*(sum(psi_e).^(-1)));
    end; %renormalization?
end;

%Return data from diagonalization for electrons
PPe=psi_e;

%------END DYNAMICS------------------------------------------------

%------OUTPUT VISUALIZATION----------------------------------------

%Wavepacket in energy space
tgrid_small=0; % time variable for energy space output
deh=0;
for tt=0:timeSTEPS-1;
t=1+tt; %*5
tgrid_small=[tgrid_small t];
slice=zeros(size(Egrid));
deh_tmp=0;

%Output population distribution as function of energy with Gaussian smear
%Electrons
for i=1:Omax-Omin+1;  
slice=slice+PPe(i,t)*exp(-((Egrid-Ee(i))/Width).^2)/(2*3.1415926*Width^2)^0.5;
end;
se=slice;

%Calculate <DeltaE>(t)
for i=1:length(Egrid);
   deh_tmp=deh_tmp+se(i)*Egrid(i)*de;
end;

s=[s slice'];

%Apply scaling factor
if t==1
  stmp=(Ee(inib))/deh_tmp; 
  ss= stmp; 
end;

deh_tmp=deh_tmp*ss;       % (average) electron hole energy at time t

deh=[deh deh_tmp];

end;%t

k=zeros(size(time));

for i=1:timeSTEPS-1;
ktmp=log((deh(2)-deh(timeSTEPS+1))/(deh(i+1)-deh(timeSTEPS+1)));    
k(i)=1/time(i)*ktmp;
end;

%Output relaxation of E(eV) to file
fid=fopen('e_t','w');
fprintf(fid,'relaxation of E(eV) \n') 
for i=1:timeSTEPS;
fprintf(fid,'%i %16.8f %16.8f \n', time(i), deh(i+1),k(i));
end;
fclose(fid);

%Ouput lowest excitation level population as a function of time to file
fid1=fopen('s1_vs_t','w');
for i=1:timeSTEPS;
fprintf(fid1,'%16.8f %16.8f %16.8f \n', time(i), PPe(1,i));
end;
fclose(fid1);

time1=[min(time) time];
%tgrid=[0 tgrid]
%mesh(tgrid_small,Egrid,s);axis([0 200 -2 1.8 -5 5]);view(2)
%mesh(tgrid_small,Egrid,s);view(2);       %population distribution

%----CHANGE HERE FOR ELECTRONS vs HOLES spatial distribution vs time---

%consider electrons spatial extent
%meshc((PPe(:,:)'*Be)');view(2);     %P_CB(z,t), only consider rho_ii(t)
%WPe=(((PPe(:,:)'*Be)'));            % WPe(z,t)

%consider holes spatial extent
%meshc((PPh(:,:)'*Bh)');view(2);     %P_VB(z,t), only consider rho_ii(t)
%WPe=(((PPh(:,:)'*Bh)'));            % WPh(z,t)


%meshc(Zz,tgrid,WPe');view(2);

%transfer=(sum(WPe(55:135,:)));           % sum over 55~135 for each t
%norm=sum(WPe);
%norm1=norm.^(-1);
%plot(transfer);
%plot(norm);
%plot(transfer.*norm);
%plot(transfer.*norm.^(-1));
%CT=[time' (transfer.*norm.^(-1))'];
%save CT CT -ASCII

%zSTEPS=600;
%Zgrid=1:zSTEPS;
%Zzgrid=max(Zz)*Zgrid/zSTEPS;           %redefine
%slice1=zeros(size(Zgrid));
%WWPe=slice1';
%ZWidth=Zz(2)-Zz(1);
%for i=1:timeSTEPS;
%    slice1=zeros(size(Zgrid));
%for j=1:zmax;
%    slice1=slice1+WPe(j,i)*exp(-((Zzgrid-Zz(j))/ZWidth).^2)/.5;
%end;
%WWPe=[WWPe slice1'];
%end;
%time1=[time(1) time];
%mesh(tgrid_small,Zzgrid,WWPe);axis([.001 400 0 30.97]);colorbar;view(2)
tgrid_small_two=zeros(1,timeSTEPS+1);
for i=1:timeSTEPS;
    tgrid_small_two(1,i+1) = time(1,i);
end;
mesh(tgrid_small,Egrid,s);axis([.001 600 -0.1 3.5]);set(gca,'FontSize',20,'FontWeight','bold');xlabel('Time (ps)','FontSize',24,'FontWeight','normal');ylabel('Energy (eV)','FontSize',24,'FontWeight','normal');view(2);caxis([0.01 5])       %population distribution
