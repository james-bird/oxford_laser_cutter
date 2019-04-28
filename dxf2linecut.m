%% ------------------------------------------------------------- %%
%          .DXF to gCode converter for the Oxford Laser           %
%                          James Bird                             %
%                    j.bird12@imperial.ac.uk                      %
%                             2015                                %
% --------------------------------------------------------------- %

% This code is designed to turn drawings created in CREO (and exported as 
% .dxf files) into cutting paths for the Oxford Laser. It uses a 
% combination of moving the bed and the galvos to maximise cutting
% speed over a larger area. The code can handle circles, arcs, 
% lines, and splines (but they must be exported as polylines in CREO). 
% The code removes the boarders from the .dxf files, and applies
% an offset. It might be necessary to apply an
% additional offset manually, if the cutting paths are in the 2nd to 4th
% quadrants, or if many cutting paths coincide with overlap regions. 
% The script is not very efficient, so feel to make changes,
% and let me know if you find any bugs!



%% Dock all figures
set(0,'DefaultFigureWindowStyle','docked');

%% Clear everything
home;
clearvars;
close all;

%% Settings/parameters

file_name='square.dxf';         %Input file name
out_name='laser_code.pgm';      %Output file names
remove_border =true;            %Remove the boarder from CREO drawing
offsett=[1,1];                  %Manually apply offset (mm) (recommended to prevent repeated cuts along the axes)
galvorange=13;                  %+/- range of the galvo (mm)
overlap=0.05;                   %Overlap of the galvo regions (mm) (recommended to account for non-vertical cutting laser)
SF=1;                           %SF of the G-code (Applied this scale-factor to g-code output - useful to massively reduce scale of part for debugging)
delay=0.1;                      %Delay applied when changing motion - 0.1s seemed to work well
bedS=40;                        %Speed of the bed
galvoS=400;                     %Speed of the galvos
nopass=50;                      %Number of double passes


%% Read in the file

[c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = f_LectDxf(file_name);  

%% Get/sort data

%Lines
minlinex=[];
minliney=[];
lines=[];
if ~cellfun('isempty',c_Line)
  lines=zeros(size(c_Line,1),4);
  for count=1:size(lines,1)
    lines(count,:)=c_Line{count,1};
  end
  if remove_border == true
    lines(1:4,:)=[];
  end
  minlinex=min(min(lines(:,[1,3])));
  minliney=min(min(lines(:,[2,4])));
end

%Polylines
polylines=[];
minpolyx=[];
minpolyy=[];
if ~cellfun('isempty',c_Poly)
  polytemp=[];
  polylines=cell(size(c_Poly,1),1);
  for count=1:size(c_Poly,1)
    polylines{count}=c_Poly{count,1};
    polytemp=[polytemp;c_Poly{count,1}];
  end
  if ~isempty(polytemp)
    minpolyx=min(polytemp(:,1));
    minpolyy=min(polytemp(:,2));
  end
end

%Circles
circles=[];
if ~cellfun('isempty',c_Cir)
  circles=zeros(size(c_Cir,1),3);
  for count=1:size(c_Cir,1)
    circles(count,:)=[c_Cir{count,1}];
  end
end

%Arcs
arcs=[];
if ~cellfun('isempty',c_Arc)
  arcs=zeros(size(c_Arc,1),5);
  for count=1:size(c_Arc,1)
    arcs(count,:)=[c_Arc{count,1}];
    if arcs(count,4)>arcs(count,5)
      arcs(count,5)=arcs(count,5)+360;
    end
  end
end

%% Find offset such that (0,0) in the bottem left of the drawing (ignoring circles and arcs)

xoff=min([minlinex,minpolyx])-offsett(1);
yoff=min([minliney,minpolyy])-offsett(2);

%Apply offset;

if ~cellfun('isempty',c_Poly)
  for count=1:size(c_Poly,1)
    ptemp=polylines{count};
    ptemp(:,1)=ptemp(:,1)-xoff;
    ptemp(:,2)=ptemp(:,2)-yoff;
    polylines{count}=ptemp;
  end
end

if ~cellfun('isempty',c_Line)
  lines(:,[1,3])=lines(:,[1,3])-xoff;
  lines(:,[2,4])=lines(:,[2,4])-yoff;
end

if ~cellfun('isempty',c_Cir)
  circles(:,1)=circles(:,1)-xoff;
  circles(:,2)=circles(:,2)-yoff;
end

if ~cellfun('isempty',c_Arc)
  arcs(:,1)=arcs(:,1)-xoff;
  arcs(:,2)=arcs(:,2)-yoff;
end

%% Break drawing down into areas

galvocount=1;
galvo.lines=[];
galvo.stuff=[];
galvo.circles=[];
galvo.arcs=[];
galvo.poly=[];
for xcount=galvorange:galvorange*2:300
  for ycount=galvorange:galvorange*2:300
    galvo(galvocount).stuff=0;    
    galvo(galvocount).xoff=xcount;
    galvo(galvocount).yoff=ycount;
    polysquare=([[xcount-(galvorange+overlap);xcount-(galvorange+overlap);xcount+(galvorange+overlap);xcount+(galvorange+overlap);xcount-(galvorange+overlap)],[ycount-(galvorange+overlap);ycount+(galvorange+overlap);ycount+(galvorange+overlap);ycount-(galvorange+overlap);ycount-(galvorange+overlap)]]);
    %Lines
    for linecount=1:size(lines,1)
      [interx,intery] = polyxpoly(lines(linecount,[1,3]),lines(linecount,[2,4]),polysquare(:,1),polysquare(:,2),'unique');
      [in,on] = inpolygon(lines(linecount,[1,3]),lines(linecount,[2,4]),polysquare(:,1),polysquare(:,2));
      insid=in;
      if size(interx,1)==2
        galvo(galvocount).lines=[galvo(galvocount).lines;[interx(1),intery(1),interx(2),intery(2)]];
      elseif sum(insid)==2
        galvo(galvocount).lines=[galvo(galvocount).lines;lines(linecount,:)];
      elseif size(interx,1)==1 && insid(1)==1
        galvo(galvocount).lines=[galvo(galvocount).lines;[lines(linecount,[1:2]),interx,intery]];
      elseif size(interx,1)==1 && insid(2)==1
        galvo(galvocount).lines=[galvo(galvocount).lines;[interx,intery,lines(linecount,[3:4])]];
      end
      if ~isempty(galvo(galvocount).lines)
        galvo(galvocount).polysquare=polysquare;
        galvo(galvocount).stuff=1;
        galvo(galvocount).lineoff=galvo(galvocount).lines;
        galvo(galvocount).lineoff(:,[1,3])=galvo(galvocount).lineoff(:,[1,3])-xcount;
        galvo(galvocount).lineoff(:,[2,4])=galvo(galvocount).lineoff(:,[2,4])-ycount;
      end
    end
    %Circles
    for circount=1:size(circles,1)
      tempcir=circles(circount,3)*exp(1i*[0:(2*pi)/1e3:2*pi]);
      cirx=circles(circount,1)+(real(tempcir));
      ciry=circles(circount,2)+(imag(tempcir));
      [centerinsid] = inpolygon(circles(circount,1),circles(circount,2),polysquare(:,1),polysquare(:,2));
      [ptsinside] = inpolygon(cirx,ciry,polysquare(:,1),polysquare(:,2));
      if sum (ptsinside)==size(cirx,2)
        galvo(galvocount).circles=[galvo(galvocount).circles;[circles(circount,[1:2]),circles(circount,1),circles(circount,2)+circles(circount,3),circles(circount,1),circles(circount,2)+circles(circount,3)]]; 
      else
        ind1=find(diff([ptsinside])==1);
        ind2=find(diff([ptsinside])==-1);
        ind1(ind1==size(cirx,2))=1;
        ind2(ind2==size(cirx,2))=1;
        for bount=1:size(ind1,2)
          [strtx,strty]=polyxpoly(cirx([ind1(bount),ind1(bount)+1]),ciry([ind1(bount),ind1(bount)+1]),polysquare(:,1),polysquare(:,2),'unique');
          [ennx,enny] = polyxpoly(cirx([ind2(bount),ind2(bount)+1]),ciry([ind2(bount),ind2(bount)+1]),polysquare(:,1),polysquare(:,2),'unique');
          galvo(galvocount).circles=[galvo(galvocount).circles;[circles(circount,1:2),strtx,strty,ennx,enny]];
        end
      end
      if ~isempty(galvo(galvocount).circles)
        galvo(galvocount).polysquare=polysquare;
        galvo(galvocount).stuff=1;
        galvo(galvocount).circoff=galvo(galvocount).circles;
        galvo(galvocount).circoff(:,[1,3,5])=galvo(galvocount).circoff(:,[1,3,5])-xcount;
        galvo(galvocount).circoff(:,[2,4,6])=galvo(galvocount).circoff(:,[2,4,6])-ycount;
      end
    end
    %Arcs
    for arccount=1:size(arcs,1)
      startt=arcs(arccount,4);
      enn=arcs(arccount,5);
      temparcs=arcs(arccount,3)*exp(1i*(pi/180)*[startt:(enn-startt)/1000:enn]);
      arcx=arcs(arccount,1)+real(temparcs);
      arcy=arcs(arccount,2)+imag(temparcs);
      [centerinsid] = inpolygon(arcs(arccount,1),arcs(arccount,2),polysquare(:,1),polysquare(:,2));
      [ptsinside] = inpolygon(arcx,arcy,polysquare(:,1),polysquare(:,2));
      if sum (ptsinside)==size(arcx,2)
        galvo(galvocount).arcs=[galvo(galvocount).arcs;[arcs(arccount,1:2),arcx(1),arcy(1),arcx(end),arcy(end)]];
      else
        ind1=find(diff([ptsinside])==1);
        ind2=find(diff([ptsinside])==-1);
        ind1(ind1==size(arcx,2))=1;
        ind2(ind2==size(arcy,2))=1;
        [stptinsid] = inpolygon(arcx(1),arcy(1),polysquare(:,1),polysquare(:,2));
        [enntinsid] = inpolygon(arcx(end),arcy(end),polysquare(:,1),polysquare(:,2));
        if stptinsid>0
          [ennx,enny] = polyxpoly(arcx([ind2(1),ind2(1)+1]),arcy([ind2(1),ind2(1)+1]),polysquare(:,1),polysquare(:,2),'unique');
          galvo(galvocount).arcs=[galvo(galvocount).arcs;[arcs(arccount,1:2),arcx(1),arcy(2),ennx,enny]];
          ind2(1)=[];
        end
        if enntinsid>0
          [strtx,strty] = polyxpoly(arcx([ind1(1),ind1(1)+1]),arcy([ind1(1),ind1(1)+1]),polysquare(:,1),polysquare(:,2),'unique');
          galvo(galvocount).arcs=[galvo(galvocount).arcs;[arcs(arccount,1:2),strtx,strty,arcx(end),arcy(end)]];
          ind1(1)=[];
        end
        for bount=1:size(ind1,2)
          [strtx,strty]=polyxpoly(arcx([ind1(bount),ind1(bount)+1]),arcy([ind1(bount),ind1(bount)+1]),polysquare(:,1),polysquare(:,2),'unique');
          [ennx,enny] = polyxpoly(arcx([ind2(bount),ind2(bount)+1]),arcy([ind2(bount),ind2(bount)+1]),polysquare(:,1),polysquare(:,2),'unique');
          galvo(galvocount).arcs=[galvo(galvocount).arcs;[arcs(arccount,1:2),strtx,strty,ennx,enny]];
        end 
      end
      if ~isempty(galvo(galvocount).arcs)
        galvo(galvocount).polysquare=polysquare;
        galvo(galvocount).stuff=1;
        galvo(galvocount).arcoff=galvo(galvocount).arcs;
        galvo(galvocount).arcoff(:,[1,3,5])=galvo(galvocount).arcoff(:,[1,3,5])-xcount;
        galvo(galvocount).arcoff(:,[2,4,6])=galvo(galvocount).arcoff(:,[2,4,6])-ycount;
      end
    end
    %Polylines
    linecount=1;
    for polycount=1:size(polylines,1)
      temppoly=[polylines{polycount}];
      [in,on] = inpolygon(temppoly(:,1),temppoly(:,2),polysquare(:,1),polysquare(:,2));
      insid=in';
      ind1=find(diff([0,insid,0])==1);
      ind2=find(diff([0,insid,0])==-1)-1;
      for bount=1:size(ind1,2)
        tempvec=[];
        if ind1(bount)>1
          if insid(ind1(bount)-1)==0
            [interx,intery] = polyxpoly(temppoly([ind1(bount),ind1(bount)-1],1),temppoly([ind1(bount),ind1(bount)-1],2),polysquare(:,1),polysquare(:,2),'unique');
            tempvec=[tempvec;[interx,intery]];
          end
        end
        tempvec=[tempvec;[temppoly(ind1(bount):ind2(bount),:)]];
       if ind2(bount)<size(temppoly,1)
          if insid(ind2(bount)+1)==0
            [interx2,intery2] = polyxpoly(temppoly([ind2(bount),ind2(bount)+1],1),temppoly([ind2(bount),ind2(bount)+1],2),polysquare(:,1),polysquare(:,2),'unique');
            tempvec=[tempvec;[interx2,intery2]];
          end
        end
        galvo(galvocount).poly{linecount}=tempvec;
        linecount=linecount+1;
      end    
    end
    if ~isempty(galvo(galvocount).poly)
      galvo(galvocount).polysquare=polysquare;
      galvo(galvocount).stuff=1;
      for dount=1:size(galvo(galvocount).poly,2)
        yemp=galvo(galvocount).poly{dount};
        yemp(:,1)=yemp(:,1)-xcount;
        yemp(:,2)=yemp(:,2)-ycount;
        galvo(galvocount).polyoff{dount}=yemp;
       end
     end
    galvocount=galvocount+1;
  end
end


%% Plot Things
%plot lines
hold on
title('Whole drawing')
xlabel('X axis (mm)');
ylabel('Y axis (mm)');
if ~cellfun('isempty',c_Line)
  for count=1:size(lines,1)
    plot(lines(count,[1,3]),lines(count,[2,4]),'b') 
  end
end
hold off

%plot polylines
if ~cellfun('isempty',c_Poly)
  hold on
  for count=1:size(polylines,1)
    templines=polylines{count};
    plot(templines(:,1),templines(:,2),'r')  
  end
  hold off
end

%plot circles
if ~cellfun('isempty',c_Cir)
  hold on
  for count=1:size(circles,1)
    templines=circles(count,3)*exp(1i*[0:0.01:2*pi]);
    plot(circles(count,1)+real(templines),circles(count,2)+imag(templines),'k')  
  end
  hold off
end

%plot arcs
if ~cellfun('isempty',c_Arc)
  hold on
  for count=1:size(arcs,1)
    startt=arcs(count,4);
    enn=arcs(count,5);
    if enn < startt
      enn=enn+360;
    end
    templines=arcs(count,3)*exp(1i*(pi/180)*[startt:0.01:enn]);
    plot(arcs(count,1)+real(templines),arcs(count,2)+imag(templines),'g')  
  end
  hold off
end
axis equal






%% Plot galvos areas

figure;
hold on
for counter=1:length(galvo)
  if galvo(counter).stuff==1
  axis equal
  xlabel('X axis (mm)');
  ylabel('Y axis (mm)');
  if ~isempty(galvo(counter).lineoff)
    for count=1:size(galvo(counter).lineoff,1)
      plot(galvo(counter).lineoff(count,[1,3])+galvo(counter).xoff,galvo(counter).lineoff(count,[2,4])+galvo(counter).yoff,'r--^') 
      plot(galvo(counter).polysquare(:,1),galvo(counter).polysquare(:,2),'-.c');
    end
  end
  
  if ~isempty(galvo(counter).circles)
    for count=1:size(galvo(counter).circoff,1)
      plot(galvo(counter).circoff(count,[1,3,5])+galvo(counter).xoff,galvo(counter).circoff(count,[2,4,6])+galvo(counter).yoff,'d-')
%       plot(galvo(counter).circles(count,[1,3,5]),galvo(counter).circles(count,[2,4,6]),'kd')
      plot(galvo(counter).polysquare(:,1),galvo(counter).polysquare(:,2),'-.g');
    end
  end
  if ~isempty(galvo(counter).arcs)
    for count=1:size(galvo(counter).arcoff,1)
%       plot(galvo(counter).arcoff(count,[1,3,5])+galvo(counter).xoff,galvo(counter).arcoff(count,[2,4,6])+galvo(counter).yoff,'d-')
      plot(galvo(counter).arcs(count,[1,3,5]),galvo(counter).arcs(count,[2,4,6]),'-d')
      plot(galvo(counter).polysquare(:,1),galvo(counter).polysquare(:,2),'-.y');
    end
  end
  if ~isempty(galvo(counter).poly)
    for count=1:length(galvo(counter).poly)
        tpoly=galvo(counter).polyoff{count};
%       plot(galvo(counter).arcoff(count,[1,3,5])+galvo(counter).xoff,galvo(counter).arcoff(count,[2,4,6])+galvo(counter).yoff,'d-')
      plot(tpoly(:,1)+galvo(counter).xoff,tpoly(:,2)+galvo(counter).yoff,'-d')
      plot(galvo(counter).polysquare(:,1),galvo(counter).polysquare(:,2),'-.k');

    end
  end
  end
end
if ~cellfun('isempty',c_Cir)
  for count=1:size(circles,1)
    templines=circles(count,3)*exp(1i*[0:0.01:2*pi]);
    plot(circles(count,1)+real(templines),circles(count,2)+imag(templines),'k')  
  end
end
if ~cellfun('isempty',c_Arc)
  for count=1:size(arcs,1)
    templines=arcs(count,3)*exp(1i*(pi/180)*[arcs(count,4):(arcs(count,5)-arcs(count,4))/100:arcs(count,5)]);
    plot(arcs(count,1)+real(templines),arcs(count,2)+imag(templines),'--g')  
  end
end
hold off





%% Create Gcode

 fid=fopen(out_name,'w');
fprintf(fid,'DVAR $DISTZ\n');
% fprintf(fid,'DVAR $MILLSTEP\n');

%Write Header
fprintf(fid,';\n; G-code to cut out lattice structure\n;\n;**************************************\n');
% fprintf(fid,'');
fprintf(fid,'G90 ; use absolute coordinate system\n');
fprintf(fid,'G71 ; units in mm\n');
fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n');
fprintf(fid,'G92 X0 Y0 Z0 ;Create new absolute datum\n');
fprintf(fid,'$DISTZ=0 \n');
fprintf(fid,'FARCALL "ATTENUATOR.PGM" s100 ; 100power\n');
fprintf(fid,'G16 GX GY Z \n');
fprintf(fid,'F%.0f \n',galvoS);

for gcount=1:size(galvo,2)
  if galvo(gcount).stuff==1;
    fprintf(fid,'DWELL %f\nG1 X%f Y%f F%.0f \nDWELL %f\nF%.0f \n',delay,galvo(gcount).xoff*SF,galvo(gcount).yoff*SF,bedS,delay,galvoS);   
    
    
    %cut out lines
    if ~isempty(galvo(gcount).lines)
    for lincount=1:size(galvo(gcount).lines,1)
      fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n');
      fprintf(fid,'G1 GX %f GY %f \n',galvo(gcount).lineoff(lincount,1)*SF,galvo(gcount).lineoff(lincount,2)*SF);
      fprintf(fid,'BEAMON ; turn beam on\n');
      fprintf(fid,'REPEAT %.0f\n',nopass);
      fprintf(fid,'G1 GX %f GY %f \n',galvo(gcount).lineoff(lincount,3)*SF,galvo(gcount).lineoff(lincount,4)*SF);
      fprintf(fid,'G1 GX %f GY %f \n',galvo(gcount).lineoff(lincount,1)*SF,galvo(gcount).lineoff(lincount,2)*SF);
      fprintf(fid,'ENDREPEAT\n');
      fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n'); 
    end
    end
    
%    cut out polylines
     if ~isempty(galvo(gcount).poly)
     for polycount=1:size(galvo(gcount).poly,2)
       temppoly=galvo(gcount).polyoff{polycount};
       fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n');
       fprintf(fid,'G1 GX %f GY %f \n',temppoly(1,1)*SF,temppoly(1,2)*SF);
       fprintf(fid,'BEAMON ; turn beam on\n');
       fprintf(fid,'REPEAT %.0f\n',nopass);
       for polcounter=2:size(temppoly,1)
         fprintf(fid,'G1 GX %f GY %f \n',temppoly(polcounter,1)*SF,temppoly(polcounter,2)*SF);   
       end
       for polcounter=size(temppoly,1)-1:-1:1
         fprintf(fid,'G1 GX %f GY %f \n',temppoly(polcounter,1)*SF,temppoly(polcounter,2)*SF);   
       end
       fprintf(fid,'ENDREPEAT\n');
      fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n'); 
     end
     end
    
     %cut out circles
     if ~isempty(galvo(gcount).circles)
       for circout=1:size(galvo(gcount).circles,1)
        cent=galvo(gcount).circoff(circout,[1,2]);
        P1=galvo(gcount).circoff(circout,[3,4]);
        P2=galvo(gcount).circoff(circout,[5,6]);
%         
          fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n');
           fprintf(fid,'G1 GX %f GY %f \n',P1(1)*SF,P1(2)*SF);
          fprintf(fid,'BEAMON ; turn beam on\n');
          fprintf(fid,'REPEAT %.0f\n',nopass);
          fprintf(fid,'G03 GX %f GY %f I %f J%f\n',P2(1)*SF,P2(2)*SF,(cent(1)-P1(1))*SF,(cent(2)-P1(2))*SF); %anti-clockwise
          fprintf(fid,'G02 GX %f GY %f I %f J%f\n',P1(1)*SF,P1(2)*SF,(cent(1)-P2(1))*SF,(cent(2)-P2(2))*SF); %clockwise
      fprintf(fid,'ENDREPEAT\n');
      fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n'); 
      end
     end
    
          %cut out arcs
     if ~isempty(galvo(gcount).arcs)
       for arccout=1:size(galvo(gcount).arcs,1)
        cent=galvo(gcount).arcoff(arccout,[1,2]);
        P1=galvo(gcount).arcoff(arccout,[3,4]);
        P2=galvo(gcount).arcoff(arccout,[5,6]);
          fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n');
           fprintf(fid,'G1 GX %f GY %f \n',P1(1)*SF,P1(2)*SF);
          fprintf(fid,'BEAMON ; turn beam on\n');
          fprintf(fid,'REPEAT %.0f\n',nopass);
          fprintf(fid,'G03 GX %f GY %f I %f J%f\n',P2(1)*SF,P2(2)*SF,(cent(1)-P1(1))*SF,(cent(2)-P1(2))*SF); %anti-clockwise
          fprintf(fid,'G02 GX %f GY %f I %f J%f\n',P1(1)*SF,P1(2)*SF,(cent(1)-P2(1))*SF,(cent(2)-P2(2))*SF); %clockwise
      fprintf(fid,'ENDREPEAT\n');
      fprintf(fid,'BEAMOFF ; Ensure that the beam is off\n'); 
      end
    end
  end
end
fprintf(fid,'DWELL %f\nG1 X%f Y%f F%.0f \nDWELL %f\nF%.0f \n',delay,0,0,bedS,delay,galvoS);
fprintf(fid,'G1 GX %f GY %f \n',0,0);
fprintf(fid,'M2 ; End the program\n');
fclose(fid);













