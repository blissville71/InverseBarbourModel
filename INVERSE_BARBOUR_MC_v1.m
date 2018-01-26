function INVERSE_BARBOUR_MC_v1(INPUT,OUTPUT,SIMNUM,MACHINE_ERROR)
%
%Code name: INVERSE_BARBOUR_MC_v1
%Authors: Michael Singer; Cristina Evans; Chris Sargeant, last updated, 2018
%Date created: 13/11/2014
%
%Inverse Barbour model to compute source water based on measured tree ring
%cellulose with error bars determined by Monte Carlo sampling of climate
%and physiological variables.
%
%Based on Barbour, M. M., J. S. Roden, G. D. Farquhar, and J. R. Ehleringer. 2004. 
%Expressing leaf water and cellulose oxygen isotope ratios as enrichment above source
%water reveals evidence of a Péclet effect. Oecologia 138:426-435, 10.2307/40005850.
%
%USAGE: INVERSE_BARBOUR_MC_v1(INPUT,OUTPUT,SIMNUM,MACHINE_ERROR)
%
%INPUT is a string containing the name of the input textfile
%
%OUTPUT is a string containing the name of the output textfile
%
%SIMNUM is number of Monte Carlo simulations to run by varying input parameters
%
%MACHINE_ERROR is decimal containing the average analytical error report from the IRMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%model parameters

%number of Monte Carlo simulations
sim_num = SIMNUM;
machine_error = MACHINE_ERROR; %enter the appropriate value here for analytical error based on machine output
confidence_lim = 1.96; %based on 95% confidence. use 2.575 for 99%; 1.96 for 95%; 1.645 for 90%

%Standard Mean Ocean Water d18O
SMOW=0.0020052;

%concentration of d18O in water (mol m-3)
concwater=55500;

%diffusivity of d18O in water (m2 s-1)
diffhwater=2.66e-9;

%diffusive fractionation of d18O through stomata (per mil)
diffstomata=28;%changed 32 to 28 following Luz 2009 study that merlivat values more accurate than cappa 2003

%diffusive fractionation through boundary layer (per mil)
diffboundary=19;%changed from 21 to 19 as above

%equilibrium fractionation between C=O and water (per mil)
equilfrac=27;
%----------------------------------------------------------------------------
%open input table

A=importdata(INPUT);
data=A.data;
header=A.textdata;
%----------------------------------------------------------------------------
%determine the accuracy of the approximation

choose = 2;
chosenaccuracy=choose+2;
for p=1:chosenaccuracy
    accuracy(p)=10^(p-1);
end
[j,~] = size(data);

%assign site number

for n=1:j
   
    %choose a starting value for the source water
    %lets make it the actual cellulose value
    
    source=data(n,1);
    
    %measured d18O in tree ring cellulose
    actualcellulose=data(n,1);
    
    %atmospheric water vapour d18O (per mil)
    atmwat_real=data(n,2);
    atmwat_low=data(n,3);
    atmwat_high=data(n,4);
    atmwat_range = (atmwat_low:0.1:atmwat_high);
    atmwat_rangey = atmwat_range;
    for boo = 1:100
        atmwat_rangey = [atmwat_rangey,atmwat_range]; %this just repeats the sequence many times to allow datasample to function better.
    end
    atmwat_range = atmwat_rangey';
    
    
    %relative humidity (%)
    relhum_real=data(n,5);
    relhum_low=data(n,6);
    relhum_high=data(n,7);
    relhum_range = (relhum_low:0.1:relhum_high);
    relhum_rangey = relhum_range;
    for boo = 1:100
        relhum_rangey = [relhum_rangey,relhum_range]; %this just repeats the sequence many times to allow datasample to function better.
    end
    relhum_range = relhum_rangey';
    
    %air temp (oC)
    airtemp_real=data(n,8);
    airtemp_low=data(n,9);
    airtemp_high=data(n,10);
    airtemp_range = (airtemp_low:0.1:airtemp_high);
    airtemp_rangey = airtemp_range;
    for boo = 1:100
        airtemp_rangey = [airtemp_rangey,airtemp_range]; %this just repeats the sequence many times to allow datasample to function better.
    end
    airtemp_range = airtemp_rangey';
    
    %leaf temperature difference (oC)
    leafdiff=1;
    
    %     %barometric pressure (kPa)
    %     barpress_real=data(n,11);
    %     barpress_low=data(n,12);
    %     barpress_high=data(n,13);
    %     barpress_range = (barpress_low:barpress_high);
    %     barpress_rangey = barpress_range;
    %     for boo = 1:100
    %         barpress_rangey = [barpress_rangey,barpress_range]; %this just repeats the sequence many times to allow datasample to function better.
    %     end
    %     barpress_range = barpress_rangey';
    
    %stomatal conductance (mol m-2 s-1)
    stomcond_real=data(n,14);
    stomcond_low=data(n,15);
    stomcond_high=data(n,16);
    stomcond_range = (stomcond_low:0.1:stomcond_high);
    stomcond_rangey = stomcond_range;
    for boo = 1:100
        stomcond_rangey = [stomcond_rangey,stomcond_range]; %this just repeats the sequence many times to allow datasample to function better.
    end
    stomcond_range = stomcond_rangey';
    
    
    
    %transpiration rate (mmol m-2 s-1)
    transp_real=data(n,17);
    transp_low=data(n,18);
    transp_high=data(n,19);
    transp_range = (transp_low:0.1:transp_high);
    transp_rangey = transp_range;
    for boo = 1:100
        transp_rangey = [transp_rangey,transp_range]; %this just repeats the sequence many times to allow datasample to function better.
    end
    transp_range = transp_rangey';
    
    %proportion of exchangeable oxygen in cellulose via Roden et al., 2000
    exchoxy=0.42;
    
    %proportion of xylem water in meristem
    propmeristem=1;
    %----------------------------------------------------------------------------
    %calculate new values under Monte Carlo scheme sampling climatic and
    %tree variables randomly for 1000 simulations in order to create error
    %estimates around our computed values of back-calculated source waters
    
    for MC = 1:1:sim_num
        
        transp = datasample(transp_range,1);
        
        
        stomcond = datasample(stomcond_range,1);
        
        
        %barpress = datasample(barpress_range,1);
        
        
        airtemp = datasample(airtemp_range,1);
        
        
        relhum = datasample(relhum_range,1);
        
        
        atmwat = datasample(atmwat_range,1);
        
        %effective path length (m) via Song et al., 2014
        transp2 = transp/1000;
        pathlength=0.0000236*transp2^(-1.2);
        path_error = 0.19*pathlength; %this error corresponds to the R2 value reported in the paper (0.81).
        path_range = pathlength-path_error:0.0001:pathlength+path_error;
        path_rangey = path_range;
        for boo = 1:100
            path_rangey = [path_rangey,path_range]; %this just repeats the sequence many times to allow datasample to function better.
        end
        path_range = path_rangey';
        pathlength = datasample(path_range,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UNFINISHED
        %boundary layer conductance (mol m-2 s-1)
        %boundcond=1;
        boundcond_range = 0.2:0.1:3;
        boundcond_rangey = boundcond_range;
        for boo = 1:100
            boundcond_rangey = [boundcond_rangey,boundcond_range]; %this just repeats the sequence many times to allow datasample to function better.
        end
        boundcond_range = boundcond_rangey';
        boundcond = datasample(boundcond_range,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %leaf temperature
        leaftemp=airtemp+leafdiff;
        
        %equilibrium vapour pressure fractionation (per mil) Majoube,1971
        equilvapfrac=((exp((1137/((airtemp+273.16)^2))-(0.4156/(airtemp+273.16))-0.0020667))-1)*1000;
        
        %total diffusive fractionation (per mil)
        totdifffrac=(((1/stomcond)*diffstomata)+((1/boundcond)*diffboundary))/((1/boundcond)+(1/stomcond));
        
        %18O/16O water vapour (per mil)
        ratiovapour=((atmwat/1000)+1)*SMOW;
        
        %18O/16O source water (per mil)
        ratiosource=((source/1000)+1)*SMOW;
        
        %Delta18O water vapour (per mil)
        Deltavapour=((ratiovapour/ratiosource)-1)*1000;
        
        %saturated vapour pressure of the air (mbar)
        satvappress=(6.13753*exp(airtemp*((18.564-(airtemp/254.4)))/(airtemp+255.57)));
        
        %actual vapour pressure of the air (mbar)
        actvappress=(relhum/100)*satvappress;
        
        %leaf internal vapour pressure (mbar)
        leafvappress=(6.13753*exp(leaftemp*((18.564-(leaftemp/254.4)))/(leaftemp+255.57)));
        
        %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
        CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
        
        %Peclet number
        peclet=((transp/1000)*pathlength)/(concwater*diffhwater);
        
        %bulk leaf water (Delta18O) (per mil)
        bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
        
        %sucrose Delta18O (per mil)
        %sucroseD=bulkleafD+equilfrac;
        
        %cellulose Delta18O (per mil)
        celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
        
        %Craig-Gordon evaporation site water (d18Oe) (per mil)
        %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %bulk leaf water (d18O) (per mil)
        %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %sucrose d18O (per mil)
        %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %cellulose d18O (per mil)
        cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
        %----------------------------------------------------------------------------
        
        for y=1:chosenaccuracy
            
            %create a loop to calculate a more accurate value
            if mod(y,2)~=0
                while cellulose*accuracy(y)>actualcellulose*accuracy(y)
                    source=source-(1/accuracy(y));
                    %18O/16O source water (per mil)
                    ratiosource=((source/1000)+1)*SMOW;
                    
                    %Delta18O water vapour (per mil)
                    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                    
                    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                    
                    %bulk leaf water (Delta18O) (per mil)
                    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                    
                    %sucrose Delta18O (per mil)
                    %sucroseD=bulkleafD+equilfrac;
                    
                    %cellulose Delta18O (per mil)
                    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                    
                    %Craig-Gordon evaporation site water (d18Oe) (per mil)
                    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %bulk leaf water (d18O) (per mil)
                    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %sucrose d18O (per mil)
                    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %cellulose d18O (per mil)
                    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                end
            else
                while cellulose*accuracy(y)<actualcellulose*accuracy(y)
                    source=source+(1/accuracy(y));
                    
                    %18O/16O source water (per mil)
                    ratiosource=((source/1000)+1)*SMOW;
                    
                    %Delta18O water vapour (per mil)
                    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                    
                    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                    
                    %bulk leaf water (Delta18O) (per mil)
                    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                    
                    %sucrose Delta18O (per mil)
                    %sucroseD=bulkleafD+equilfrac;
                    
                    %cellulose Delta18O (per mil)
                    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                    
                    %Craig-Gordon evaporation site water (d18Oe) (per mil)
                    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %bulk leaf water (d18O) (per mil)
                    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %sucrose d18O (per mil)
                    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %cellulose d18O (per mil)
                    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                end
            end
        end
        
        %----------------------------------------------------------------------------
        %store values to plot in a graph
        
        sourcestor(MC,n)=source;
        
    end
    std_source(n) = std(sourcestor(:,n));
    [~,~] = size(sourcestor);
    se_source(n) = 2*confidence_lim*std_source(n)./sqrt(sim_num); % confidence value defined at top of code
    %maxsource(n) = max(sourcestor(:,n));
    %minsource(n) = min(sourcestor(:,n));
    %----------------------------------------------------------------------------
    %redo the calculations for the measured (mean) values for climatic and
    %physiological variables
    atmwat=atmwat_real;
    transp=transp_real;
    airtemp=airtemp_real;
    stomcond=stomcond_real;
    relhum=relhum_real;
    %barpress=barpress_real;
    
    %effective path length (m) via Song et al., 2014
    pathlength=0.0000236*transp^(-1.2);
    
    %leaf temperature
    leaftemp=airtemp+leafdiff;
    
    %equilibrium vapour pressure fractionation (per mil)
    equilvapfrac=((exp((1137/((airtemp+273.16)^2))-(0.4156/(airtemp+273.16))-0.0020667))-1)*1000;
    
    %total diffusive fractionation (per mil)
    totdifffrac=(((1/stomcond)*diffstomata)+((1/boundcond)*diffboundary))/((1/boundcond)+(1/stomcond));
    
    %18O/16O water vapour (per mil)
    ratiovapour=((atmwat/1000)+1)*SMOW;
    
    %18O/16O source water (per mil)
    ratiosource=((source/1000)+1)*SMOW;
    
    %Delta18O water vapour (per mil)
    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
    
    %saturated vapour pressure of the air (mbar)
    satvappress=(6.13753*exp(airtemp*((18.564-(airtemp/254.4)))/(airtemp+255.57)));
    
    %actual vapour pressure of the air (mbar)
    actvappress=(relhum/100)*satvappress;
    
    %leaf internal vapour pressure (mbar)
    leafvappress=(6.13753*exp(leaftemp*((18.564-(leaftemp/254.4)))/(leaftemp+255.57)));
    
    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
    
    %Peclet number
    peclet=((transp/1000)*pathlength)/(concwater*diffhwater);
    
    %bulk leaf water (Delta18O) (per mil)
    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
    
    %sucrose Delta18O (per mil)
    %sucroseD=bulkleafD+equilfrac;
    
    %cellulose Delta18O (per mil)
    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
    
    %Craig-Gordon evaporation site water (d18Oe) (per mil)
    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %bulk leaf water (d18O) (per mil)
    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %sucrose d18O (per mil)
    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %cellulose d18O (per mil)
    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
    %----------------------------------------------------------------------------
    
    for y=1:chosenaccuracy
        
        %create a loop to calculate a more accurate value
        if mod(y,2)~=0
            while cellulose*accuracy(y)>actualcellulose*accuracy(y)
                source=source-(1/accuracy(y));
                %18O/16O source water (per mil)
                ratiosource=((source/1000)+1)*SMOW;
                
                %Delta18O water vapour (per mil)
                Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                
                %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                
                %bulk leaf water (Delta18O) (per mil)
                bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                
                %sucrose Delta18O (per mil)
                %sucroseD=bulkleafD+equilfrac;
                
                %cellulose Delta18O (per mil)
                celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                
                %Craig-Gordon evaporation site water (d18Oe) (per mil)
                %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %bulk leaf water (d18O) (per mil)
                %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %sucrose d18O (per mil)
                %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %cellulose d18O (per mil)
                cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
            end
        else
            while cellulose*accuracy(y)<actualcellulose*accuracy(y)
                source=source+(1/accuracy(y));
                
                %18O/16O source water (per mil)
                ratiosource=((source/1000)+1)*SMOW;
                
                %Delta18O water vapour (per mil)
                Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                
                %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                
                %bulk leaf water (Delta18O) (per mil)
                bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                
                %sucrose Delta18O (per mil)
                %sucroseD=bulkleafD+equilfrac;
                
                %cellulose Delta18O (per mil)
                celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                
                %Craig-Gordon evaporation site water (d18Oe) (per mil)
                %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %bulk leaf water (d18O) (per mil)
                %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %sucrose d18O (per mil)
                %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %cellulose d18O (per mil)
                cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
            end
        end
    end
    sourcestore(n)=source;
    cellulosestore(n)=cellulose;
    actualcellulosestore(n)=actualcellulose;
    
end
%----------------------------------------------------------------------------
%save the output
%labels = {header{2:j+1,1}}';
%labels = header(2:j+1,1)';
%labels = char(header{2:j+1,1});
source = sourcestore;
% maxsource = maxsource;
% minsource = minsource;
upper_source = source+se_source+machine_error; %added in the average analytical error here defined above.
lower_source = source-se_source-machine_error;
Out=[sourcestore;upper_source;lower_source];
%Out=[labels;source;maxsource;minsource];

%plot the results
figure
x = 1:j;
X=[x,fliplr(x)];
Y=[upper_source,fliplr(lower_source)];
fill(X,Y,'g')
hold on
h = plot(source,'r');
%Fig = get(h, 'Children');
set(h, 'LineWidth', 0.5)
ylabel('Modeled Source Water d18O (per mil)');
xlabel('Samples plotted as a time series');
hold off
M=Out';
%open the output file
inverseoutput=fopen(OUTPUT,'w');

fprintf(inverseoutput,'%11s\t %25s\t %21s\t %21s\t %22s\t %28s\n','Sample_Code','Modeled_Source_Water_d18O','Source_Water_SE+_d18O', 'Source_Water_SE-_d18O', 'Modeled Tree Ring d18O', 'Measured Tree d18O');
%fprintf(inverseoutput,'%6s\t %6.2f\t %6.2f\t %6.2f\t\n',Out);
%fprintf(inverseoutput,'%s\n',header{2:j+1,1});
for row = 1:j
    %fprintf(inverseoutput,'%s\t %6.2f\t %6.2f\t %6.2f\n',header{row+1,1}, M(row,:));
    fprintf(inverseoutput,'%s\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\n',header{row+1,1}, M(row,:), cellulosestore(row), actualcellulosestore(row));
end
%fprintf(inverseoutput,'%s\t %6.2f\t 6.2f\t 6.2f\n', header{2:j+1,1}, M(1:j,:));
%fprintf(inverseoutput,'%s\n',header{end});
%fprintf(inverseoutput, [repmat('\t %6.2f\t',1,size(M,2)-1) '%6.2f\n'], M');
%fprintf(inverseoutput,'%6.2f\t %6.2f\t %6.2f\t\n',Out);
fclose(inverseoutput);

%disp('min')