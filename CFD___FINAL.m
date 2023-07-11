%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project Computational Fluid Dynamics
% Profossor : DR.Naderan 
% Student : Vahid Eftekhari Khorasani
% Student Number: 401126104
% Amirkabir University of Technology
% Wedge flow for compressible flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc ;
clear ;  

%%%%%%%%%%%%%%%%%%
% number of Mesh %
%%%%%%%%%%%%%%%%%%
degree = 1;  % 1 for 15 degre and 2 for 35 degree

if degree ~= 3    
    N_m=5;
else
    N_m=1;    
end        
Fact = zeros(4,N_m);
U_Mean = zeros(4,N_m);
H = zeros(1,N_m); 
Error_1=zeros(4,N_m,N_m);


%%%%%%%%%%%%%%%%%%%%%%%%%
% FREESTREAM CONDITIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%

T_inf        = 288;                            % Kelvin
P_inf        = 101325;                         % Freestream pressure [Pa]
M_inf        = 2.5;                            % Freestream Mach number
R_air        = 287;                            % Specific gas constant for air at standard conditions [J*kg^-1*K^-1]
Theta        = 0;                              % deg -- Angle of attack (AOA)
rho_inf      = P_inf/(R_air*T_inf);            % kg/m^3

a_inf        = SpeedOfSound(P_inf,rho_inf);    % m/s
u_inf        = M_inf*a_inf*cosd(Theta);        % m/s
v_inf        = M_inf*a_inf*sind(Theta);        % m/s
vel_inf      = [u_inf v_inf];

% Freestream Primitive State Vector (PSV)
V_free        = [rho_inf u_inf v_inf P_inf];

diary('output.txt')
fprintf('Freestream Conditions:\n');
fprintf('Mach:         %10.2f\n',M_inf);
fprintf('Flow AOA:     %10.2f deg\n',Theta);
fprintf('u Velocity:   %10.2f m/s\n',u_inf);
fprintf('v Velocity:   %10.2f m/s\n',v_inf);
fprintf('Pressure:     %10.2e Pa\n',P_inf);
fprintf('Temperature:  %10.2f K\n',T_inf);
fprintf('Density:      %10.2f kg/m^3\n\n',rho_inf);

%%%%%%%%%%%%%%%%%%%%%%%%
%  ITERATION VARIABLES %
%%%%%%%%%%%%%%%%%%%%%%%%

iterations   = 1000;            % Number of iterations
                               
timestep     = 1e-5;            % Timestep for global timestepping
CFL          = 0.5;             % Courant number CFL

t_stage      = 4;               % m-stage time stepping
                                % e.g. 1 for Euler step
                                % 4 for 4th-order RK


% Output variables
fprintf('Iteration Variables:\n');
fprintf('Iterations:   %5d\n',iterations);
fprintf('CFL:          %5.2f\n',CFL);
fprintf('M-Stage:      %5d\n',t_stage);


for aa=1:N_m    

residual={zeros(1,1)}  ; V={zeros(1,1)} ; U={zeros(1,1)};
len_x_e=zeros(1,1)     ; len_y_e=zeros(1,1);
len_x_n=zeros(1,1)     ; len_y_n=zeros(1,1);
len_x_w=zeros(1,1)     ; len_y_w=zeros(1,1);
len_x_s=zeros(1,1)     ; len_y_s=zeros(1,1);
x_mid=zeros(1,1)       ; y_mid=zeros(1,1);
volume=zeros(1,1)      ; mov=zeros(1,1);
sE=zeros(1,1)     ; sN=zeros(1,1) ; sW=zeros(1,1) ; sS=zeros(1,1);
nE={zeros(1,1)}     ; nN={zeros(1,1)} ; nW={zeros(1,1)} ; nS={zeros(1,1)};

%%%%%%%%%%%%%%%%%%%%%%%%
% GRID SIZE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

if degree==1
    if aa==1
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge15_1.x');
    elseif aa==2
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge15_2.x');
    elseif aa==3
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge15_3.x');
    elseif aa==4
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge15_4.x');
    else
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge15_5.x');
    end

elseif degree==2
    if aa==1
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge35_1.x');
    elseif aa==2
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge35_2.x');
    elseif aa==3
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge35_3.x');
    elseif aa==4
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge35_4.x');
    else
    fid = fopen('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD - 2\project_2\Geometry\Wedge35_5.x');
    end
else
    disp('Geometry Not founded')
end
    
    if fid >= 1
        % Read in file headers
        zones = fscanf(fid, '%d', 1);
        % Code only handles 1 zone
        % Therefore, check for number of zones
        if (zones == 1)
            % Read in number of i,j,k points
            npi = fscanf(fid, '%d', 1);
            npj = fscanf(fid, '%d', 1);
            npk = fscanf(fid, '%d', 1);
       
            % Retrieve i,j,k coordinates
            x = fscanf(fid, '%f', [npi,npj]);
            y = fscanf(fid, '%f', [npi,npj]);
            z = fscanf(fid, '%f', [npi,npj]);
            disp('Grid read successfully');
            rho_y=zeros(npi-1,N_m);
            U_y=zeros(npi-1,N_m);
            V_y=zeros(npi-1,N_m);
            P_y=zeros(npi-1,N_m);
            X_y=zeros(npi-1,N_m);
            Fac_y=zeros(npi-1);
            Fac_x=zeros(npi-1);
          
        end
        fclose(fid);
    end

    %%%%%%%%%%%%%%%%
    % GRID METRICS %
    %%%%%%%%%%%%%%%%

z_width = [0,0,1];                   % Unit vector in z-dir
nc_i = npi-1;                         % Number of cells (pts-1) in i dir
nc_j = npj-1;                         % Number of cells (pts-1) in j dir

for i = 1:nc_i
    for j = 1:nc_j
        % Assemble the face lengths
        len_x_e(i,j) = x(i+1,j+1)-x(i+1,j);
        len_y_e(i,j) = y(i+1,j+1)-y(i+1,j);
        
        len_x_n(i,j) = x(i,j+1)-x(i+1,j+1);
        len_y_n(i,j) = y(i,j+1)-y(i+1,j+1);

        len_x_w(i,j) = x(i,j)-x(i,j+1);
        len_y_w(i,j) = y(i,j)-y(i,j+1);

        len_x_s(i,j) = x(i+1,j)-x(i,j);
        len_y_s(i,j) = y(i+1,j)-y(i,j);
        
        % Compute midpoint of cell (for plotting)
        x_mid(i,j) = (x(i,j) + x(i+1,j))/2;
        y_mid(i,j) = (y(i,j) + y(i,j+1))/2;
        
        % Compute volume of cell using A.BxC (volume of parallelepiped)
        volume(i,j) = abs(dot(z_width,cross(-1*[len_x_s(i,j),len_y_s(i,j),0],...
                      [len_x_e(i,j),len_y_e(i,j),0])));
        
        % Compute area of cell
        sE(i,j) = sqrt((len_x_e(i,j))^2 + (len_y_e(i,j))^2);
        sN(i,j) = sqrt((len_x_n(i,j))^2 + (len_y_n(i,j))^2);
        sW(i,j) = sqrt((len_x_w(i,j))^2 + (len_y_w(i,j))^2);
        sS(i,j) = sqrt((len_x_s(i,j))^2 + (len_y_s(i,j))^2);
        
        % Compute outward normal of faces (return 3 component vector)
        temp_nE = cross([len_x_e(i,j),len_y_e(i,j), 0]/sE(i,j), z_width);
        temp_nN = cross([len_x_n(i,j),len_y_n(i,j), 0]/sN(i,j), z_width);
        temp_nW = cross([len_x_w(i,j),len_y_w(i,j), 0]/sW(i,j), z_width);
        temp_nS = cross([len_x_s(i,j),len_y_s(i,j), 0]/sS(i,j), z_width);
        
        % Truncate normal vector to 2 components
        nE{i,j} = [temp_nE(1) temp_nE(2)];
        nN{i,j} = [temp_nN(1) temp_nN(2)];
        nW{i,j} = [temp_nW(1) temp_nW(2)];
        nS{i,j} = [temp_nS(1) temp_nS(2)];
        
        
    end 
end


%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

resid_ii      = 0;               % Iterative residual
resid_a      = 0;               % Step 0 residual
start_it1   = 0;               % Used for multiple runs
end_iter     = 0;               % Used for multiple runs
residualReduced = 0;               % If divergence detected
fm = 1;                         % Used for capturing movies

% Combine normals and areas into big cell array which will be passed
% to the function which computes the residual
normals = {nE nN nW nS};
areas   = {sE sN sW sS};

% Initalize variables which will allow for visualization
% i.e. Plot Contours
m_density = zeros([nc_i,nc_j]);
m_uvel    = zeros([nc_i,nc_j]);
m_vvel    = zeros([nc_i,nc_j]);
m_pres    = zeros([nc_i,nc_j]);

% Loop through all cells and init PSV to freestream conditions
% Convert PSV to conservative state vector
% Init residual to 0
for i = 1:nc_i
    for j = 1:nc_j
        V{i,j}     = V_free;
        U{i,j}     = consvV_U(V{i,j});
        residual{i,j} = [0 0 0 0];
    end
end


    %%%%%%%%%%%%%
    % MAIN LOOP %
    %%%%%%%%%%%%%
start_it1 = start_it1 + 1;        % Start at iteration 1

% Main loop in time
for iter = start_it1:(end_iter + iterations)
    % Time variable used to measure time/iteration
    ti1 = cputime;
    
    % Initialize iteration residual to 0
    resid_ii = 0;
    
    % Save CSV from this timestep (to be used in m-stage)
    U0 = U;
    
    % M-stage timestepping scheme
    for m = 1:t_stage
        % Calculate residual using function calcResid
        % Passes PSV, normals, areas cell array,
        %        freestream PSV, and nci and ncj
        residual = calcResid(V, V_free, normals, areas, nc_i, nc_j);
        
        % Loop through all cells to update solution
        for i = 1:nc_i
            for j = 1:nc_j
                
                
                    vel_inf = [V{i,j}(2) V{i,j}(3)];
                    cell_a = SpeedOfSound(V{i,j}(4),V{i,j}(1));
                    dt(1) = CFL * sE(i,j)/(abs(vel_inf(1)*nE{i,j}(1) +...
                            vel_inf(2)*nE{i,j}(2))+cell_a);
                    dt(2) = CFL * sN(i,j)/(abs(vel_inf(1)*nN{i,j}(1) +...
                            vel_inf(2)*nN{i,j}(2))+cell_a);
                    dt(3) = CFL * sW(i,j)/(abs(vel_inf(1)*nW{i,j}(1) +...
                            vel_inf(2)*nW{i,j}(2))+cell_a);
                    dt(4) = CFL * sS(i,j)/(abs(vel_inf(1)*nS{i,j}(1) +...
                            vel_inf(2)*nS{i,j}(2))+cell_a);
                    timestep = min(dt);
            
                
                % Update solution using the saved CSV
                % Multiply by 'alpha' constant 
                U{i,j} = U0{i,j} - 1/(t_stage-(m-1))*timestep/volume(i,j)*residual{i,j};
                
                % Update cell PSV
                V{i,j} = consvU_V(U{i,j});
            
                % Update contour arrays used for plotting
                m_density(i,j) = V{i,j}(1);
                m_uvel(i,j)    = V{i,j}(2);
                m_vvel(i,j)    = V{i,j}(3);
                m_pres(i,j)    = V{i,j}(4);

                % Assemble first part of L2 norm for residual
                resid_ii = resid_ii + (residual{i,j}.^2);
            end
        end
    end   
    
    % Assemble second part of L2 norm for residual
    resid_ii = (resid_ii).^.5/(nc_i*nc_j);
    
    % Assign normalization value in first interation
    if iter == 1
        resid_a = resid_ii;
         
    end
    
    % Detects extreme divergence (at the point of no return)
    % and shuts down simulation
    if isnan(resid_ii/resid_a)
        disp('Solution corrupt.');
        break;
    end
    
    % Detects divergence happening in x-mom resid and cuts CFL in half
        ee=(resid_ii(2)/resid_a(2));
        if ee >= (1e1)
            if residualReduced == 0
                CFL = CFL/2; 
                notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
                disp(notice);
                residualReduced = residualReduced + 1;
            end
        end

    % Computes time/iteration
    ti2 = cputime-ti1;
    Error_1(:,iter,aa)=[resid_ii(1)/resid_a(1);resid_ii(2)/resid_a(2);resid_ii(3)/resid_a(3);resid_ii(4)/resid_a(4)];
       
  if(max(abs(Error_1(:,iter,aa))))<=1e-4 
      break;
  end
  if(isreal(m_density))==0 ||  max(abs(Error_1(:,iter)))>200
      disp('The solution is divergent')
      break;
  end
end

start_it1 = iter;
end_iter = iter;

%% Mach matrix calculation
 h=mean(mean(len_y_e));

    %%%%%%%%%%%%%%%
    % POSTPROCESS %
    %%%%%%%%%%%%%%%

    % Density (kg/m^3) contour
    figure(1);
%     subplot(2,3,ii);
    surf(x_mid',y_mid',m_density');
    shading interp
    colorbar
    title([' Density (kg/m^3) (M_{\infty} = ' , num2str(M_inf),')' ])
    xlabel('x')
    ylabel('y')
    view(2)

    % Velocity X (m/s) contour
    figure(2);
%     subplot(2,3,ii);
    surf(x_mid',y_mid',m_uvel');
    shading interp
    colorbar
    title([' Velocity X (m/s) (M_{\infty} = ' , num2str(M_inf),')' ])
    xlabel('x')
    ylabel('y')
    view(2)

    % Velocity Y (m/s) contour
    figure(3);
%     subplot(2,3,ii);
    surf(x_mid',y_mid',m_vvel');
    shading interp
    colorbar 
    title([' Velocity Y (m/s) (M_{\infty} = ' , num2str(M_inf),')' ])
    xlabel('x')
    ylabel('y')
    view(2)
 
    % Pressure (Pa)
    figure(4);
%     subplot(2,3,ii);
    surf(x_mid',y_mid',m_pres');
    shading interp
    colorbar 
    title([' Pressure (pa) (M_{\infty} = ' , num2str(M_inf),')' ])
    xlabel('x')
    ylabel('y')
    view(2)

    


    % Error 
    ite=1:iter;
    figure(10);
    subplot(2,3,aa)
    plot(ite,Error_1(1,1:iter,aa),'r')
    hold on
    plot(ite,Error_1(2,1:iter,aa),'b')
    plot(ite,Error_1(3,1:iter,aa),'g')
    plot(ite,Error_1(4,1:iter,aa),'k')
    legend('cont. resid','x-mom resid','y-mom resid','energy resid')
    xlabel('Iterations')
    ylabel('Residual Error')
    title(sprintf('Residual Error '))

    yy=abs(y-0.5);
    xx=abs(x-0.5);
   
    for i=1:nc_i
       Fac_y(i)=find(yy(i,:)==min(yy(i,:)),1);
       Fac_x(i)=find(xx(:,i)==min(xx(:,i)),1);
    end
    for i=1:nc_i
        U_y(i,aa)=m_uvel(i,Fac_y(i));
        V_y(i,aa)=m_vvel(i,Fac_y(i));
        P_y(i,aa)=m_pres(i,Fac_y(i));
        rho_y(i,aa)=m_density(i,Fac_y(i));
        X_y(i)= x(i,Fac_y(i));
    end
    

    figure(5);
    subplot(2,3,aa);
    plot(X_y',U_y(:,aa)');
    grid on
    xlabel('x')
    ylabel('y')
    title(sprintf('U(x) in y=0.5  '))

    figure(6);
    subplot(2,3,aa);
    plot(X_y',P_y(:,aa)');
    grid on
    xlabel('x')
    ylabel('y')
    title(sprintf('Pressure(x) in y=0.5  '))
    

    figure(7);
    subplot(2,3,aa);
    plot(X_y',V_y(:,aa)');
    grid on
    title(sprintf('V(x) in y=0.5   '))

    figure(8);
    subplot(2,3,aa);
    plot(X_y',rho_y(:,aa)');
    grid on
    xlabel('x')
    ylabel('y')
    title(sprintf('Density(x) in y=0.5  '))
    
    H(aa)=h;
    Fact(:,aa)=[ nc_j/2 ;nc_j*(5/8); nc_j*(3/4);nc_j*(7/8)];
    Fact=floor(Fact);
for i=1:4
    if degree==1
        U_Mean(i,aa)=m_uvel(nc_i*7/8,Fact(i,aa));
    else
        U_Mean(i,aa)=m_uvel(floor((npi+1)/2),Fact(i,aa));
    end

end


end

if N_m == 5 
Er=Error(U_Mean,N_m);
delta_y=Error_Slope( Er,H,N_m );

%%%%%Draw the relative error value in terms of H
COL=['s','*','d','o'];
for i=1:4
figure(12)
subplot(1,2,1)
loglog(H(1:end-1),abs(Er(i,:)),sprintf(COL(i)),'linewidth',2);
title('Solution convergence ');
ylabel('successive error');
xlabel('h');
hold on
grid on

%%%%Draw the error slope in terms of h
subplot(1,2,2)
loglog(H(1:end-2),abs(delta_y(i,:)),'linewidth',2);
title('Error slope');
ylabel('Slope');
xlabel('h');
hold on
grid on
end
legend('y=1.5','y=1.875','y=2.25','y=2.625')
end
