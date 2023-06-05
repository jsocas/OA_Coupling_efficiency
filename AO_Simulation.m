clear all;
%%
%Fiber Coupling Efficiency

% General Params
% atm realisations
iterations = 1;

nPx = 360;
D = 0.1;
Robs = 0.34;
nLenslet = 12; % DM 11x11
% ITEFI r0 = 1cm @ 850 nm -> 0.8cm @ 550nm
r0 = 0.008;
% AO freq
sampling_time = 1/2000;
% Wavelengths
lambda_comms = photometry.FSOC_1550;
lambda_beacon = photometry.R;
link_distance = 2e3;

% Fibre coupling params
w0=10e-6;
% D/f optimum from coupling_Ruilier
D_f = 0.195;


% Closed loop params don't need to be defined each time
loopGain = 0.5;
loopIter = 2000;
loopIter_noCOrrection = 1000;

% DM doesn't need to be defined each iteration
comm_signal = source('wavelength',lambda_comms,'height',link_distance); %GEO 35876000
% Satellite beacon to measure
beacon = source('wavelength',lambda_beacon,'height',link_distance); %Beacon at 850nm
tel = telescope(D,...
    'resolution',nPx,...
     'obstructionRatio', Robs,...
    'samplingTime',sampling_time);
wfs = shackHartmann(nLenslet,nPx,0.8);
beacon = beacon.*tel*wfs;
wfs.INIT

bif = influenceFunction('monotonic',50/100);
nActuator = nLenslet + 1;
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);

beacon=beacon.*tel;
dmCalib = calibration(dm,wfs,beacon,beacon.wavelength/40); %we can use a calibration source at 1550nm?
dmCalib.nThresholded = 10;
commandMatrix = dmCalib.M;

% Params for fiber coupling computation
lambda = comm_signal.wavelength;
% focal optima D/f = 0.22
f= tel.D/D_f;
% Crear mesh conveniente con el tamaño tel.D para las
% variables X e Y para el back propagated mode (fórmula Chen).
comm_signal = comm_signal.*tel;
EA = comm_signal.amplitude.*exp(1i.*comm_signal.phase);
sizeEA = size(EA,1);
X = (1:1:sizeEA) - 0.5*sizeEA;
X = tel.D*X/sizeEA;
Y = X';

FCE_atm = zeros(iterations,loopIter);

for i=1:iterations

    clear ngs atm tel wfs

    % Random seed
    s = randi([1 1000],1);
    s = RandStream('mt19937ar', 'seed', s);

    % Satellite downlink to correct
    comm_signal = source('wavelength',lambda_comms,'height',link_distance); %GEO 35876000

    % Satellite beacon to measure
    beacon = source('wavelength',lambda_beacon,'height',link_distance); %Beacon at 850nm

    % ATMOSPHERE
    atm = atmosphere(photometry.V,r0,30,...
        'altitude',[0,0.5,1,2]*1e3,...
        'fractionnalR0',[0.25,0.25,0.25,0.25],...
        'windSpeed',[5,10,20,5],...
        'windDirection',[0,pi/4,pi,pi],...
        'randStream',s);

%     atm = atmosphere(photometry.V,r0,30,...
%         'altitude',[0,4,10]*1e3,...
%         'fractionnalR0',[0.7,0.25,0.05],...
%         'windSpeed',[5,10,20],...
%         'windDirection',[0,pi/4,pi]);

%     atm = ogsAtmosphere(10,1,'semilla',s);

    % TELESCOPE
    tel = telescope(D,...
        'resolution',nPx,...
        'obstructionRatio', Robs,...
        'samplingTime',sampling_time);




    % WFS
    wfs = shackHartmann(nLenslet,nPx,0.8);
    %
    beacon = beacon.*tel*wfs;
    wfs.INIT

    %
    tel = tel+atm;

    % Resetting the DM command
    dm.coefs = 0;

    figure(11)
%     subplot(1,2,1); imagesc(comm_signal.meanRmPhase);
%     subplot(1,2,2); imagesc(beacon.meanRmPhase);
%     drawnow;

    for k=1:loopIter

        sprintf('ATM iter: %d; AO loop iter: %d', i, k)

        if k<=loopIter_noCOrrection
            +tel;
            comm_signal=comm_signal.*tel;
            beacon = beacon.*tel;
            dm.coefs = 0;
            comm_signal=comm_signal*dm;
            beacon = beacon*dm*wfs;
        else
            % Propagation throught the atmosphere to the telescope, +tel means that
            % all the layers move of one step based on the sampling time and the
            % wind vectors of the layers
            +tel;
            beacon = beacon.*tel*dm*wfs;
            comm_signal=comm_signal.*tel*dm;
%             % Saving the turbulence aberrated phase
%             turbPhase = beacon.meanRmPhase;
% %             % Variance of the atmospheric wavefront
% %             total(k) = var(ngs);
%             % Propagation to the WFS
%             comm_signal=comm_signal*dm;
%             beacon = beacon*dm*wfs;
%             % Variance of the residual wavefront
%             residue(k) = var(ngs);
            % Computing the DM residual coefficients
            residualDmCoefs = commandMatrix*wfs.slopes;
            % Integrating the DM coefficients
            dm.coefs = dm.coefs - loopGain*residualDmCoefs;
            % Display of turbulence and residual phase
%             gcf;
%             subplot(1,2,1); imagesc(comm_signal.meanRmPhase);
%             subplot(1,2,2); imagesc(beacon.meanRmPhase);
% % % set(h,'Cdata',[comm_signal.meanRmPhase,beacon.meanRmPhase])
%             drawnow;
        end

        %% Chen
        %     % Aplicación de la fórmula encontrada en el artículo de Chen.
        %     % http://dx.doi.org/10.1364/AO.54.008722
        %     %
        %     % Sin atmósfera
        %     tel = tel - atm;
        %     ngs = ngs.*tel*wfs;
        %
        %     % Onda propagada
        %     % amplitude = ngs.amplitude;
        %     % phase = ngs.phase;
        %     EA = ngs.amplitude.*exp(1i.*ngs.phase); % Dimensiones 250 x 250
        %     % % EA está "sampleado" en las dimensiones de la apertura del telescopio de
        %     % % tel.D de tamaño.

        EA_atm = comm_signal.amplitude.*exp(1i.*comm_signal.meanRmPhase);



        %     [FCE(i),~] = eta_Chen(EA,w0,lambda,f,tel.D);
        FCE_atm(i,k) = eta_Chen(EA_atm,w0,lambda,f,tel.D);


    end





end

save('May23_ITEFI_FCE_WFS633_DM1550_1atm_100AOoff_2900AOon_12x12_sampling1_2000_D25_4cm', 'FCE_atm')

figure()
% plot([1:iterations],FCE);
% hold on;
plot([1:loopIter],FCE_atm);
% hold off;
% legend('w/o atmosphere', 'w/ atmosphere');
ylabel('Coupling Efficiency \eta');
xlabel('AO Loop Iterations');


%%
%
% figure()
% plot(tel.D./f,FCE_atm./max(FCE_atm));
% hold on;
% plot(tel.D./f,etaF_atm_check./max(etaF_atm_check),'+');
% hold off;
% % %% focal fija, diferentes atmósferas
% % focal_opt = tel.D/ra_max_atm;
% % iterations = 10;
% % etaF_EA1 = zeros(1,iterations);
% % for i = 1:iterations
% %     ngs = ngs.*+tel*camH;
% %     EA1 = ngs.amplitude.*exp(1i.*ngs.phase);
% %     [etaF_EA1(i),~] = eta_Chen(EA1,w0,lambda,focal_opt,tel.D);
% % end
% % plot(1:iterations,etaF_EA1);
% %
% % k=2*pi/lambda;
%
% %% Diferentes "pasos de atmósfera" con diferentes focales
% Natm = 10;
% etaF_EA1 = zeros(Natm,size(f,2));
% for j = 1:Natm
%     for k = 1:Natm
%         ngs = ngs.*+tel*wfs;
%     end
%     EA_atm = ngs.amplitude.*exp(1i.*ngs.phase);
%     for i = 1:size(f,2)
%
%         [etaF_EA1(j,i),~] = eta_Chen(EA_atm,w0,lambda,f(i),tel.D);
%     end
% end
% %%
% figure()
% hold on;
%
% for i =1:Natm
%     txt = num2str(i);
%     plot(tel.D./f,etaF_EA1(i,:),'DisplayName',txt);
% end
% hold off;
% legend show
%
%
% % check if wave changes
% before.amplitude = ngs.amplitude;
% before.phase = ngs.phase;
% ngs = ngs.*+tel*wfs;
% after.amplitude = ngs.amplitude;
% after.phase = ngs.phase;
% substraction = after.phase - before.phase;
%
% %% cambio manual de fase para propagar beam en otra direccion - CHECK
% beam0.phase = zeros(size(ngs.phase));
% beam0.amplitude = ngs.amplitude;
% fopt = tel.D/ra_max;
% wk = 2*pi/lambda/1e6;
% ivector = 0:wk/100:wk;
% eta_EA_beam0 = zeros(size(ivector));
% k=1;
% for index=ivector
%     ph.max = index;
%     ph.min = -ph.max;
%     phaseVar = linspace(ph.min,ph.max,size(beam0.phase,1));
%     beam0.phase = repmat(phaseVar,size(ngs.phase,1),1); % Todos los elementos de una misma columna tienen la misma fase (dirección y) que va cambiando linealmente con el cambio de columna.
%     EA_beam0 = beam0.amplitude.*exp(1i.*beam0.phase);
%
%     eta_EA_beam0(k) = eta_Chen(EA_beam0,w0,lambda,fopt,tel.D);
%     k=k+1;
% end
% figure()
% plot(ivector,eta_EA_beam0);
%
% %%
