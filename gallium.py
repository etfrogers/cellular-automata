from automata.rule import Rule
from automata.automata import Automata
from math import pi, exp
from dlmread import readnfile
import numpy as np
import matplotlib.pyplot as plt
#from scipy import interpolate
#import cmath

#translation of etr_gaus_rect_borde_unifor_deca_rais_9

lattice=0.55*1e-9;   abp2=0.5; absorb=0.5;   KB=1.38078*1e-23; h=6.6260693*1e-34; c=3*1e8;   #surface melting thickness denote the initial temperature
fluence=20*1e1;    lam = 775e-9;     f0=c/lam;    w=2*pi*f0; # make sure the smae wavelengh in function_ab      #fluence J/m^2  &&&&&&&&&&&&&&&&
timestep=1*1e-15;  tao=100*timestep;     Tao=tao/timestep;   #time step,    #pulse width tao

taoend=100*1e-7;  tao32=1*1e-12;  tao2121=1*1e-9;  # after Taoend, then all the meta will decay to ground ,
tao21=1000*1e-8;  Tao21=10*1e-9;   Taoend=taoend/timestep; #decay2121 isthe middle process decay , tao21 and Tao21 is the last decay process &&& Taoend >=ontime
tao1212=1*1e-9;
decay32=1-exp(-timestep/tao32);   decay2121=1-exp(-timestep/tao2121);    decay21=1-exp(-timestep/tao21);  Decay21=1-exp(-timestep/Tao21);
Decay213=0.01;    N3decay=24;
#Raise12=1-exp(-timestep/tao1212);
Raise12=0.4*1e-4;
raise12=Raise12;
absorption_style='small';# 'small' or 'biger'
pulse_wave='gauss_N';# 'gauss_Y'  or 'gauss_N'
decay_style='border3';#  'borderY'  or 'borderN' 'bordern' last one relati to his neibour
time_domain='gaussTY';#  'gaussTY'  or  'gaussTN'
raise_style='uniformY';# 'uniformY' or 'uniformN'
reflectivity_style='averageN';# 'averageY' or 'averageN'  #averageN is simple mode just add the reflectiveity layer let it thickness is uniform
DetN32='considerN';# 'considerY or 'considerN'; that is to say decay32 if it is dependent on the neighbour

m=100;   n=101;   n1=4;   M=m/3;   excited=0.85;  meta=0.35;  ground=0; #M is the waist of the pulse

#[ mg, ng] = meshgrid(1:n,1:m);
neighbourhood = [[1, 1, 1],
                 [1, 0, 1],
                 [1, 1, 1]]

detN32=5;    # onset from excited to meta with Neibour <detN32
detN21=6;    # onset from meta to ground with Neibour <detN21 with properbility decay21 in the whole melting process(that is 1,2,3,4)
detN12=5;    ## onset from ground to meta with Neibour >detN12 definitely
N1212=1;    #from ground to meta with propability when N2121<Neibour <detN21(that's between 2-4)

detN13=4;    ## onset from ground to excited with Neibour >detN13
N21=7;  # end from meta to ground when with chnage probability Decay21 when Neibour <N21

if (pulse_wave == 'gauss_Y'):
    gao=1
elif (pulse_wave == 'gauss_N'):
    gao=0
else:
    raise Exception
    
if (time_domain == 'gaussTY'):
    tim=1 
    Ontime=2*Tao
    power=0.9394*fluence/tao
    
elif (time_domain == 'gaussTN'):
    tim=0  
    Ontime=Tao
    power=1*fluence/tao
else:
    raise Exception
if (raise_style == 'uniformY'):
    rai=0
elif (raise_style == 'uniformN'):
    rai=1 #tim=0(1) denotes rectangular pulse(gauss pulse); gao=1,0 is gauss shape or#rectangular shap;#rai=0(1) denotes the rasi12 is constant or exponentially decrease;
else:
    raise Exception

N0=power/h/f0*timestep*(lattice)**2 #general photonic number per timestep at a certain row i


nrload = readnfile('n_poly.txt','\t');krload = readnfile('k_poly.txt','\t');   # Load Ga refractive indices and select at relevant wavelength
nSol = np.interp(lam, nrload[:,0],nrload[:,1]) + 1j*np.interp(lam, krload[:,0],krload[:,1])
nrload2 = readnfile('n_liq.txt','\t');krload2 = readnfile('k_liq.txt','\t');
nLiq = np.interp(lam, nrload2[:,0],nrload2[:,1]) + 1j*np.interp(lam, krload2[:,0],krload2[:,1])
es = nSol*nSol;    el =nLiq*nLiq;   ks=nSol*2*pi*f0/c;  kl=nLiq*2*pi*f0/c; Ls=1/np.imag(ks);   LL=1/np.imag(kl);  #pennetration depth come from the wave vector
abp1=2*lattice/Ls;# cell absorption rate for a photon
#
if (absorption_style == 'small'):
    abp=abp1
elif (absorption_style == 'biger'):
    abp=abp2
else:
    raise Exception


aut = Automata([m,n], [ground, meta, excited])
aut.doPause = False
#
#figure(1)
#plotbutton=uicontrol('style','pushbutton','string','Run','fontsize',12, 'position',[100,400,50,20], 'callback', 'run=1;'); #define the run,stop,quit and number button
#erasebutton=uicontrol('style','pushbutton','string','Stop','fontsize',12,'position',[200,400,50,20],'callback','freeze=1;');
#quitbutton=uicontrol('style','pushbutton','string','Quit','fontsize',12,'position',[300,400,50,20],'callback','stop=1;close;');
#number1 = uicontrol('style','text', 'string','1', 'fontsize',12,'position',[20,400,50,20]);##[75,368,50,20]);
#
#cells= zeros(m,n);    
aut.layout[0:m,0:n1] = meta
aut.layout[0,0] = excited
#aut.layout[1:m,n1+1:n] = ground;  # initial state for the matrix and draw

aut.show('Initial State')

#from meta to ground
aut.addRule(Rule(sFrom=meta,sTo=ground,sNb=[meta,excited],N=[0, detN21],prob=decay2121))
aut.addRule(Rule(sFrom=ground, sTo=meta, sNb=[meta,excited], N=[detN12,aut.maxNbours])) #always happens if enough neighbours
aut.addRule(Rule(sFrom=ground, sTo=meta, sNb=[meta,excited], N=[N1212,detN12], prob=raise12)) #stocahstically happens if only some neighbours
                 
NSteps = 200

time = timestep*np.array(range(0,NSteps))
LaserIntensity = N0*np.exp(-2.7726*tim*(time/tao-1)**2)
LaserIntensity[time>Ontime] = 0

#imagesc(cells); xlabel('Cells','fontsize',13),ylabel('Cells','fontsize',13);set(gca,'FontSize',12);
stop = False; run = 0; freeze = 0; k=1; 
run = 1;
stepNumber = -1; # increment at the start so move to 0 during first iteration
while (not stop):########### main loop begin
    stepNumber+=1
    print(stepNumber)
    aut.evolve()
    aut.show(str(stepNumber))
    stop=stepNumber>=NSteps-1 # stops AFTER this run, so -1
#    if (run==1) 
    #labelling
#        stepnumber = 1 + str2num(get(number1,'string'));# photonic excited at a certain matrix
#        set(number1,'string',num2str(stepnumber)); 
         #a=stepnumber-1;
#        
    #Laser excitation 
     
#        
#        
#        ################### below from ground photoinduced to exicted
#        m0=1:m;
#        if (stepNumber == 0):
#            Ni=LaserIntensity(a).*exp(-gao*(m0-m/2).^2./M.^2).* Ab;      
#        else
#            Ni=LaserIntensity(a).*exp(-gao*(m0-m/2).^2./M.^2)* absorb ;     
#        end
#        
#        proo=rand(size(m0));  #####
#        Ni(Ni<proo) = 0;
#            
#        Nz=fix(Ni);
#        Nf=Ni-Nz;##specific photonic number per timestep at a certain row i
#            
##         if(Nf(m0)>0) #ths code is almost "if true" as (x-fix(x) >0) is almost always true
##             proo=rand();
##         end    
#        proo=rand(size(m0));  #####
#        Nz(Nf>proo) = Nz(Nf>proo) + 1;
#        
#        for n0=1:n-2 ## guarrentee the last colume has no excited state
#            randno = rand(size(m0));
#            prob_abs = 1-(1-abp).^Nz; #compound probability for Nz photons that at least one is absorbed
#            excited_inds = (cells(:,n0)==ground).' & (prob_abs>randno);
#            cells(excited_inds, n0) = excited;
#            Nz(excited_inds) = Nz(excited_inds)-1;
#        end
#        
#        #####################################################################
#        ###########################excited above
#        
#        
#        
#        
#        
#        if (a<=Taoend)#############
#            
#            
#            # original code uses circular boundary conditions in m (top and
#            # botom) and crops field in n. This seems sensible.
#            ##################  excited begins to decay to meta  (they can be relative to neighbour)
#            
#            pro=rand(size(cells));
#            if strcmp(DetN32, 'considerY')
#                wrapped_cells = [cells(m,:);cells;cells(1,:)]; #add wrapping (last line above first and first line below last)
#                number = filter2(neighbourhood, (wrapped_cells == meta | wrapped_cells==excited));
#                number = number(2:end-1,:); # cut off wrapped lines;
#                
#                switched_cells = (cells==excited & number<detN32 & pro<decay32);####that's
#                
#            elseif strcmp(DetN32, 'considerN')
#                switched_cells = (cells==excited & pro<decay32);####that's
#            end
#            cells(switched_cells) = meta;
#                
#            ########### excited decay to meta above
#            
    #Rule(sFrom=meta,sTo=ground,sNb=[meta,excited],N=detN21,prob=decay2121)

    #Rule(sFrom=ground, sTo=meta, sNb=[meta,excited], N=detN12) #always happens if enough neighbours
    #Rule(sFrom=ground, sTo=meta, sNb=[meta,excited], N=[N1212,detN12], prob=raise12) #stocahstically happens if only some neighbours
#            
#            
#            
#            
#            
#            
#            if (a<=Ontime)############ below is to get how many meta and get the T to decide the next Ab
#                
#                if strcmp(reflectivity_style, 'averageY')   #########
#                    
#                    for mm=1:m 
#                        count00(mm)=0;
#                        for nn=n1+1:n 
#                            if (cells(mm,nn)==excited||cells(mm,nn)==meta)
#                                count00(mm)=count00(mm)+1 ;
#                            end
#                        end
#                        
#                        tru=ismember(excited,cells(mm,:));
#                        if (tru==1)
#                            xee(mm)=max(find(cells(mm,:) == excited));
#                        else
#                            xee(mm)=0;
#                        end
#                        xmm(mm)=max(find(cells(mm,:) == meta));
#                        ctt(mm)=max(xee(mm),xmm(mm))-n1;
#                        dbb(mm)=ctt(mm)*lattice;
#                        if (dbb(mm)==0)
#                            ebb(mm)=es;
#                        else
#                            ebb(mm)=(count00(mm)*el+(ctt(mm)-count00(mm))*es)/ctt(mm) ;
#                        end
#                        
#                        Ab(mm)=function_AbC(dbb(mm),ebb(mm))  ;
#                        
#                    end####R T Ang Absorption
#                    
#                end    ########### strcmp(reflectivity_style, 'averageY')
#                
#                
#                
#                
#                if strcmp(reflectivity_style, 'averageN')
#                    
#                    excited_cells = (cells == excited) | (cells == meta);
#                    count00 = sum(excited_cells, 2);
#                    Da=count00*lattice;
#                    #function_Ab is not vectorised yet.
#                    Ab = function_Ab(Da, lam, nSol, nLiq).';
#                                    end   ########### strcmp(reflectivity_style, 'averageN')
#                
#                
#            end  ####  <= ontime
#            
#            ##########################record meta number below
#            excited_cells_n1 = (cells(:, n1+1:end) == excited) | (cells(:, n1+1:end) == meta);
#            count0 = sum(excited_cells_n1(:));
#            
#            
#            
#        end  ####yu a< taoend yizhi  Taoend
#        #################################
#        
#        
#        
#        #####perform decay caculation after Tao21
#        if (a>Taoend) 
#            #currently never runs, as Taoend = 1e10 (timesteps)
#            #this code not tested.
#            ##################################
#            
#            #this already happens above. Should it really happen again?
#            pro = rand(size(cells));
#            cells(cells == excited & pro<decay32) = meta;
#            
#            
#            if strcmp(decay_style, 'borderY') ### border decay
#                ##########################record meta number
#                excited_cells_n1 = (cells(:, n1+1:end) == excited) | (cells(:, n1+1:end) == meta);
#                count2 = sum(excited_cells_n1(:));
#               
#                Ratio=count2/count0;
#                lastdecay21=Decay21*Ratio;
#                ############
#                
#                pro=rand(size(cells));
#           
#                wrapped_cells = [cells(m,:);cells;cells(1,:)]; #add wrapping (last line above first and first line below last)
#                number = filter2(neighbourhood, (wrapped_cells == meta | wrapped_cells==excited));
#                number = number(2:end-1,:); # cut off wrapped lines;
#                
#                switched_cells = (ng>=n+1 & ng<=n-1 & cells==meta & ...
#                                  number<N21 & pro<lastdecay21);####that's
#                
#                cells(switched_cells) = ground;
#                
#                
#            end  ############ border decay
#            
#            
#            
#            
#            if strcmp(decay_style, 'borderN') ### uniform decay
#                
#                ############decay to ground uniform
#                
#                for m0=1:m 
#                    for n0=n1+1:n-1 
#                        if(cells(m0,n0)==meta)
#                            pro=rand();
#                            if(pro< decay21)
#                                cells(m0,n0)=ground;
#                            else
#                                cells(m0,n0)=meta;
#                            end
#                        end
#                    end
#                end
#            end #####################decay uniform
#            
#            
#            
#            
#            
#            
#            if strcmp(decay_style, 'border3') ### border3 decay
#                ##########################record meta number
#                
#                
#                pro=rand(size(cells));
#           
#                wrapped_cells = [cells(m-1:m,:);cells;cells(1:2,:)]; #add wrapping (last line above first and first line below last)
#                number = filter2(neighbourhood5, (wrapped_cells == meta | wrapped_cells==excited));
#                number = number(3:end-2,:); # cut off wrapped lines;
#                
#                lastdecay21=Decay213*exp(-number/N3decay);
#                switched_cells = (cells==meta & pro<lastdecay21);####that's
#                
#                cells(switched_cells) = ground;
#                            
#                cells(switched_cells) = ground;
#                
#                
#                
#                
#            end  ############ border3 decay
#            
#            
#            
#            
#            
#            
#        end ######decay process >Taoend
#        
#        
#        
#        
#        ############################################# below is caculated different raise12
#        #raise12=Raise12*exp(-2*rai*Count/Ls)
#        raise12=Raise12*exp(-rai*a/1000);
#        
#        figure(1)
#        #subplot(2,1,1)
#        imagesc(cells)
##         xx00=pcolor(cells);
##         set(xx00, 'LineWidth',0.01); ##'1' is not a valid value. Use one of these values: '-' | '--' | ':' | '-.' | 'none'.   'LineStyle','-'
#        xlabel('Cells','fontsize',13),ylabel('Cells','fontsize',13);
#        set(gca,'FontSize',12);
#        drawnow
#        
#        
#        ####################below caculated the motel number within 10 rows
#        
#        Ninter=10;countb=0;countc=0;countd=0;counte=0;countf=0;
#        
#        #################################### Bbbbbb
#        for i=1:m 
#            for j=n1+1:n1+Ninter 
#                if (cells(i,j)==excited||cells(i,j)==meta)
#                    countb=countb+1;
#                end
#            end
#        end
#        totalb=Ninter*m;
#        eb=(countb*el+(totalb-countb)*es)/totalb ;  #composited
#        db=Ninter*lattice;
#        
#        ##########################################   Ccccccc
#        for i=1:m 
#            for j=n1+Ninter+1:n1+2*Ninter 
#                if (cells(i,j)==excited||cells(i,j)==meta)
#                    countc=countc+1;
#                end
#            end
#        end
#        totalc=Ninter*m;
#        ec=(countc*el+(totalc-countc)*es)/totalc ;  #composited
#        dc=Ninter*lattice;
#        
#        
#        ##########################################Ddddddddddd
#        for i=1:m 
#            for j=n1+2*Ninter+1:n1+3*Ninter 
#                if (cells(i,j)==excited||cells(i,j)==meta)
#                    countd=countd+1;
#                end
#            end
#        end
#        totald=Ninter*m;
#        ed=(countd*el+(totald-countd)*es)/totald ;  #composited
#        dd=Ninter*lattice;
#        
#        ##########################################  EEEEEeeeeeeeeeeeee
#        for i=1:m 
#            for j=n1+3*Ninter+1:n1+4*Ninter 
#                if (cells(i,j)==excited||cells(i,j)==meta)
#                    counte=counte+1;
#                end
#            end
#        end
#        totale=Ninter*m;
#        ee=(counte*el+(totale-counte)*es)/totale ;  #composited
#        de=Ninter*lattice;
#        
#        ##########################################  FFFFFFF ffffffffffffff
#        cells_f = cells(:, n1+4*Ninter+1:n);
#        countf = sum(sum(cells_f == excited | cells_f == meta));
#        totalf = numel(cells_f);
#        ef=(countf*el+(totalf-countf)*es)/totalf ;  #composited
#        df=(n-(n1+4*Ninter))*lattice;
#        
#        
#        #ec=ea   #composited
#        xx(k)=a;#Coun(k)=count;zz(k)=ct;
#        
#        ybr(k)=real(eb); ybi(k)=imag(eb); ddb(k)=db;
#        ycr(k)=real(ec); yci(k)=imag(ec); ddc(k)=dc;
#        ydr(k)=real(ed); ydi(k)=imag(ed); ddd(k)=dd;
#        yer(k)=real(ee); yei(k)=imag(ee); dde(k)=de;
#        yfr(k)=real(ef); yfi(k)=imag(ef); ddf(k)=df;
#        
#        k=k+1;
#        #
#        # figure(2)
#        #      subplot(2,1,2)
#        #     plot(xx,Coun,'r.','markersize',10);grid;   xlabel('time step'),ylabel('number');grid;  legend('(metallic phase number)')
#        # plot(xx,yyr,'b.',xx,Coun,'r.','markersize',10);grid;   xlabel('time step'),ylabel('magnitude');grid;  legend('real(Ec)','(intermediate number)')
#        #    axis([0 500 -60 60]);
#        
#        
#        Step=xx'; #calculation step number
#        TT=xx'*timestep; #real time
#        Db=ddb'; Ebr=ybr';Ebi=ybi';# The effective permittivity of the composited
#        Dc=ddc'; Ecr=ycr';Eci=yci';# The effective permittivity of the composited
#        Dd=ddd'; Edr=ydr';Edi=ydi';# The effective permittivity of the composited
#        De=dde'; Eer=yer';Eei=yei';# The effective permittivity of the composited
#        Df=ddf'; Efr=yfr';Efi=yfi';# The effective permittivity of the composited
#        
#        
#        qE=[Step,TT,Db,Ebr,Ebi,Dc,Ecr,Eci,Dd,Edr,Edi,De,Eer,Eei,Df,Efr,Efi];
#        
#        #save e:\zlw\gallium1.txt qE -ASCII#####
#        #save \\filestore.soton.ac.uk\users\Lz1r15\mydesktop\PABLO\gallium.txt qE -ASCII#
#        #save \\filestore.soton.ac.uk\users\Lz1r15\mydesktop\ZLW\??2016 UK\cellular automata\Gallium pahse Trasition\gallium.txt qE -ASCII#
##         save C:\users\Ed\Desktop\gallium1.txt qE -ASCII#
#        # if(stepnumber==round(Ontime))
#        # pause(1);
#        
#        #     if(a==Tao-1)
#        #       freeze=1; end
#               
#    else
#        pause(0.25);
#    end
#    
#    if (freeze==1)
#        run = 0;
#        freeze = 0;
#    end
#    # if rem(a,10)==0
#    
#    # end
#    
plt.figure(2)
plt.plot(time, LaserIntensity)
plt.show()

