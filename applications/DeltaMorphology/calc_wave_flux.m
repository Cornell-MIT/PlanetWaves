function Qsmax = calc_wave_flux(E_pdf,sang)
    % Calculates the maximum wave transport to the left or right of a shoreline by waves
    % can be compared with riverine input to get R = Q_river / Qsmax where Qsmax is this function

    % based on Nienhuis+2015
    % E_pdf is wave climate energy pdf defined over range 0:359 degrees (wave) for direction waves coming from
    % for example, 
    %   if wave = 0 then waves traveling East->West, if wave = 90 then waves traveling from
    %   North -> South
    % sang is the angle of the normal vector CCW from East of the regional shoreline without a delta
    % for example,
    % if sang = 0 then shoreline goes from North -> South with liquid on the East side
    % if sang = 90 then shoreline goes from West -> East with liquid on top


    make_main_plot = 1;
    make_subplot = 0;


    % direction waves coming from
    wave = 0:359; % degrees CCW from E

    if make_main_plot 
        figure
        polarplot(deg2rad(wave),E_pdf)
        title('direction waves coming from')
    end

    thetas = -90:90;          % all possible delta flank orientation (CCW) theta = 0 corresponds to
                              % the same orientation as the regional shoreline. positive theta 
                              % corresponds to the left flank of the delta (with theta = 90 meaning 
                              % a left delta flank that juts out perpendicular to the coast)
                              % and negative theta corresponds to the right flank of the delta
                              % (with theta = 90 meaning a right delta flank that juts out 
                              % perpendicular to the coast). The left and right flanks are defined
                              % looking towards the liquid.

   Qsmax = zeros(size(sang));

   for c = 1:numel(sang) % orientation of non-deltaic shoreline

    regional_shoreline = sang(c);

    
    phi0 = wrapTo180(wave - regional_shoreline); % Calculate phi0, the wave approach directions relative to the current
                                                      % shoreline orientation. phi0 = 0 means waves are coming straight on
                                                      % towards the coast. phi0 = 90 means waves are propagating to the
                                                      % right as you look offshore. phi0 = -90 means waves are propagating 
                                                      % to the left as you look offshore. Compare with Jaap's Fig. 1A,B.

    % Zero out portions of wave climate not seen for a given regional 
    % shoreline orientation. Should be zero for winds coming from land 
    % (zero fetch, so no waves).
    E_pdf_relative = E_pdf;
    E_pdf_relative(phi0<-90) = 0;
    E_pdf_relative(phi0>90) = 0;


    Qsnet = zeros(size(thetas));
    for i = 1:numel(thetas)  % all possible delta morphs 
    
        theta = thetas(i); % the delta shoreline orientation relative to the local shoreline. 
                           % Values between 0 and -90 are the right flank of the delta 
                           % (looking offshore), and values between 0 and 90 are the left 
                           % flank (looking offshore).
       

        % longshore transport
        LST = CERC(1, 1, wrapTo180(phi0 - theta));  
        LST(isnan(LST)) = 0;
        Qsnet(i) = sum(E_pdf_relative.*LST); % Jaap supplement Equation 4

    end

    Qsmaxright = max(Qsnet(thetas<=0)); % most positive right-flank flux
    thetamaxR = thetas(Qsnet == Qsmaxright);
    Qsminleft = min(Qsnet(thetas>=0)); % most negative left-flank flux
    thetaminL = thetas(Qsnet == Qsminleft);

    Qsmax(c) = Qsmaxright - Qsminleft;

    if make_subplot
        figure;
        tiledlayout("horizontal")
        
        nexttile
        [phi0sort, I] = sort(phi0);
        plot(phi0sort,E_pdf_relative(I),'-k') % Jaap Fig. 1B
        hold on;
        xline(wave_mean - sang(c),'-r','wave mean - regional shoreline')
        set(gca,'xdir','reverse','xlim',[-180 180])
        xlabel('Waves wrt Regional, \phi_0 = wave - regional (\circ)')
        ylabel('Wave energy PDF')
        title(['Regional shoreline angle: ' num2str(regional_shoreline) '\circ'])
    
        nexttile
        plot(thetas,CERC(1,1,thetas),'-k') % Fig. 1C
        xline(0)
        yline(0)
        grid on
        set(gca,'xdir','reverse','xlim',[-90 90])
        xlabel('Delta angle \theta (\circ)')
        ylabel('Q_s')
        title(['wave mean  - regional shore: ' num2str(wave_mean - sang(c))])
        
        nexttile
        plot(thetas,Qsnet,'-k') % Fig. 1D
        hold on
        plot(thetaminL,Qsminleft,'ob')
        plot(thetamaxR,Qsmaxright,'or')
        xline(0)
        yline(0)
        xline(wave_mean - sang(c),'-r','wave mean - regional shoreline')
        grid on
        set(gca,'xdir','reverse','xlim',[-90 90])
        xlabel('Deltaic shoreline angle, \theta (\circ)')
        ylabel('Net Q_s ')
        title(['Q_{s,max} = ' num2str(Qsmaxright) ' - ' num2str(Qsminleft) ' = ' num2str(Qsmax(c))])
    end

   end

   if make_main_plot 
       if numel(sang) == 1
           fprintf('Qs,max = %0.2f',Qsmax(1))
       else
           figure;
           plot(sang,Qsmax,'.k')
           xlabel('regional shoreline orientation')
           ylabel('Q_{s,max}')
       end
   end

end