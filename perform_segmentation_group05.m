function res=perform_segmentation(folder,segfolder)
% res=perform_segmentation(folder,segfolder)
% perform very simple segmentaiton using matlab's built in functions
% segfolder is the subfolder within folder to which we write the
% segmentation files, folder should contain jpg images to be segmented
 
imext='jpg';
segext='png';
 
D=dir(fullfile(folder,sprintf('*.%s',imext)));
if (~exist(fullfile(folder,segfolder),'dir'))
    mkdir(folder,segfolder);
end
 
for i=1:length(D)
    fullFileName = fullfile(folder,D(i).name);
    [curfolder, curfilebase, curfileext] = fileparts(fullFileName);
    segFileName = fullfile(folder,segfolder, sprintf('%s.%s',curfilebase, segext));
    timeFileName = fullfile(folder,segfolder, sprintf('%s.%s',curfilebase, 'txt'));
    % we just overwrite files for now
    % Check if files exist.
    %if exist(segFileName, 'file')
    %    errorMessage = sprintf('Error: %s exists in given folder.', segFileName);
    %    disp(errorMessage);
    %    continue;
    %end
    org=imread(fullFileName); % original image
    tic; % start timer
    seg=perform_single_segmentation(org); % this one uses grayscale image, so
    timetaken=toc; % stop timer
    
    imwrite(seg,segFileName,'png','BitDepth',1); % write result in segfolder
    fid=fopen(timeFileName,'w'); % open file for writing time
    fprintf(fid,'%.5e',timetaken); % write time taken
    fclose(fid);
end
 
end
 
function seg=perform_single_segmentation(I)
% performs segmentation on gray scale image I
% read help text for imsegfmm
 INEW=I;
 
Final=rgb2gray(I);
%Sizelar
sizeX=size(Final,2);
sizeY=size(Final,1);
 
%r-g-b densities
r = I(:,:,1);
g = I(:,:,2);
b = I(:,:,3);
cont=r==0; 
r(cont)=1;
cont=g==0; 
g(cont)=1;
cont=b==0; 
b(cont)=1;
%im to double
r=double(r);
g=double(g);
b=double(b);
 
%Density Oranlar?
rg=log(r./g);rb=log(r./b);
br=log(b./r);bg=log(b./g);
gb=log(g./b);gr=log(g./r);
 
 
Border=10;
BorderX=1;
BorderY=1;
Border2X=50;
Border2Y=50;
 
 
%ÜST
%Max ve min Diff bul
rg1Ust=rg(1:Border,BorderX:Border2X); rg2Ust=rg(2:Border+1,BorderX:Border2X); diffrgUst=rg1Ust-rg2Ust; maxdiffrgUst=max(diffrgUst(:)); mindiffrgUst=min(diffrgUst(:)); minrgUst=min(rg1Ust(:)); maxrgUst=max(rg1Ust(:)); medianrgUst=median(rg1Ust(:)); mediandiffrgUst=median(diffrgUst(:));
rb1Ust=rb(1:Border,BorderX:Border2X); rb2Ust=rb(2:Border+1,BorderX:Border2X); diffrbUst=rb1Ust-rb2Ust; maxdiffrbUst=max(diffrbUst(:)); mindiffrbUst=min(diffrbUst(:)); minrbUst=min(rb1Ust(:)); maxrbUst=max(rb1Ust(:)); medianrbUst=median(rb1Ust(:)); mediandiffrbUst=median(diffrbUst(:));
 
br1Ust=br(1:Border,BorderX:Border2X); br2Ust=br(2:Border+1,BorderX:Border2X); diffbrUst=br1Ust-br2Ust; maxdiffbrUst=max(diffbrUst(:)); mindiffbrUst=min(diffbrUst(:)); minbrUst=min(br1Ust(:)); maxbrUst=max(br1Ust(:)); medianbrUst=median(br1Ust(:)); mediandiffbrUst=median(diffbrUst(:));
bg1Ust=bg(1:Border,BorderX:Border2X); bg2Ust=bg(2:Border+1,BorderX:Border2X); diffbgUst=bg1Ust-bg2Ust; maxdiffbgUst=max(diffbgUst(:)); mindiffbgUst=min(diffbgUst(:)); minbgUst=min(bg1Ust(:)); maxbgUst=max(bg1Ust(:)); medianbgUst=median(bg1Ust(:)); mediandiffbgUst=median(diffbgUst(:));
 
gb1Ust=gb(1:Border,BorderX:Border2X); gb2Ust=gb(2:Border+1,BorderX:Border2X); diffgbUst=gb1Ust-gb2Ust; maxdiffgbUst=max(diffgbUst(:)); mindiffgbUst=min(diffgbUst(:)); mingbUst=min(gb1Ust(:)); maxgbUst=max(gb1Ust(:)); mediangbUst=median(gb1Ust(:)); mediandiffgbUst=median(diffgbUst(:));
gr1Ust=gr(1:Border,BorderX:Border2X); gr2Ust=gr(2:Border+1,BorderX:Border2X); diffgrUst=gr1Ust-gr2Ust; maxdiffgrUst=max(diffgrUst(:)); mindiffgrUst=min(diffgrUst(:)); mingrUst=min(gr1Ust(:)); maxgrUst=max(gr1Ust(:)); mediangrUst=median(gr1Ust(:)); mediandiffgrUst=median(diffgrUst(:));
 
%SOL
%Max ve min Diff bul
rg1Sol=rg(BorderY:Border2Y,1:Border); rg2Sol=rg(BorderY:Border2Y,2:Border+1); diffrgSol=rg1Sol-rg2Sol; maxdiffrgSol=max(diffrgSol(:)); mindiffrgSol=min(diffrgSol(:)); minrgSol=min(rg1Sol(:)); maxrgSol=max(rg1Sol(:)); medianrgSol=median(rg1Sol(:)); mediandiffrgSol=median(diffrgSol(:));
rb1Sol=rb(BorderY:Border2Y,1:Border); rb2Sol=rb(BorderY:Border2Y,2:Border+1); diffrbSol=rb1Sol-rb2Sol; maxdiffrbSol=max(diffrbSol(:)); mindiffrbSol=min(diffrbSol(:)); minrbSol=min(rb1Sol(:)); maxrbSol=max(rb1Sol(:)); medianrbSol=median(rb1Sol(:)); mediandiffrbSol=median(diffrbSol(:));
 
br1Sol=br(BorderY:Border2Y,1:Border); br2Sol=br(BorderY:Border2Y,2:Border+1); diffbrSol=br1Sol-br2Sol; maxdiffbrSol=max(diffbrSol(:)); mindiffbrSol=min(diffbrSol(:)); minbrSol=min(br1Sol(:)); maxbrSol=max(br1Sol(:)); medianbrSol=median(br1Sol(:)); mediandiffbrSol=median(diffbrSol(:));
bg1Sol=bg(BorderY:Border2Y,1:Border); bg2Sol=bg(BorderY:Border2Y,2:Border+1); diffbgSol=bg1Sol-bg2Sol; maxdiffbgSol=max(diffbgSol(:)); mindiffbgSol=min(diffbgSol(:)); minbgSol=min(bg1Sol(:)); maxbgSol=max(bg1Sol(:)); medianbgSol=median(bg1Sol(:)); mediandiffbgSol=median(diffbgSol(:));
 
gb1Sol=gb(BorderY:Border2Y,1:Border); gb2Sol=gb(BorderY:Border2Y,2:Border+1); diffgbSol=gb1Sol-gb2Sol; maxdiffgbSol=max(diffgbSol(:)); mindiffgbSol=min(diffgbSol(:)); mingbSol=min(gb1Sol(:)); maxgbSol=max(gb1Sol(:)); mediangbSol=median(gb1Sol(:)); mediandiffgbSol=median(diffgbSol(:));
gr1Sol=gr(BorderY:Border2Y,1:Border); gr2Sol=gr(BorderY:Border2Y,2:Border+1); diffgrSol=gr1Sol-gr2Sol; maxdiffgrSol=max(diffgrSol(:)); mindiffgrSol=min(diffgrSol(:)); mingrSol=min(gr1Sol(:)); maxgrSol=max(gr1Sol(:)); mediangrSol=median(gr1Sol(:)); mediandiffgrSol=median(diffgrSol(:)); 
 
%ALT
%Max ve min Diff bul
rg1Alt=rg(sizeY-Border-1:sizeY-1,BorderX:Border2X); rg2Alt=rg(sizeY-Border:sizeY,BorderX:Border2X); diffrgAlt=rg1Alt-rg2Alt; maxdiffrgAlt=max(diffrgAlt(:)); mindiffrgAlt=min(diffrgAlt(:)); minrgAlt=min(rg1Alt(:)); maxrgAlt=max(rg1Alt(:)); medianrgAlt=median(rg1Alt(:)); mediandiffrgAlt=median(diffrgAlt(:));
rb1Alt=rb(sizeY-Border-1:sizeY-1,BorderX:Border2X); rb2Alt=rb(sizeY-Border:sizeY,BorderX:Border2X); diffrbAlt=rb1Alt-rb2Alt; maxdiffrbAlt=max(diffrbAlt(:)); mindiffrbAlt=min(diffrbAlt(:)); minrbAlt=min(rb1Alt(:)); maxrbAlt=max(rb1Alt(:)); medianrbAlt=median(rb1Alt(:)); mediandiffrbAlt=median(diffrbAlt(:));
 
br1Alt=br(sizeY-Border-1:sizeY-1,BorderX:Border2X); br2Alt=br(sizeY-Border:sizeY,BorderX:Border2X); diffbrAlt=br1Alt-br2Alt; maxdiffbrAlt=max(diffbrAlt(:)); mindiffbrAlt=min(diffbrAlt(:)); minbrAlt=min(br1Alt(:));  maxbrAlt=max(br1Alt(:)); medianbrAlt=median(br1Alt(:)); mediandiffbrAlt=median(diffbrAlt(:));
bg1Alt=bg(sizeY-Border-1:sizeY-1,BorderX:Border2X); bg2Alt=bg(sizeY-Border:sizeY,BorderX:Border2X); diffbgAlt=bg1Alt-bg2Alt; maxdiffbgAlt=max(diffbgAlt(:)); mindiffbgAlt=min(diffbgAlt(:)); minbgAlt=min(bg1Alt(:));  maxbgAlt=max(bg1Alt(:)); medianbgAlt=median(bg1Alt(:)); mediandiffbgAlt=median(diffbgAlt(:));
 
gb1Alt=gb(sizeY-Border-1:sizeY-1,BorderX:Border2X); gb2Alt=gb(sizeY-Border:sizeY,BorderX:Border2X); diffgbAlt=gb1Alt-gb2Alt; maxdiffgbAlt=max(diffgbAlt(:)); mindiffgbAlt=min(diffgbAlt(:)); mingbAlt=min(gb1Alt(:)); maxgbAlt=max(gb1Alt(:)); mediangbAlt=median(gb1Alt(:)); mediandiffgbAlt=median(diffgbAlt(:));
gr1Alt=gr(sizeY-Border-1:sizeY-1,BorderX:Border2X); gr2Alt=gr(sizeY-Border:sizeY,BorderX:Border2X); diffgrAlt=gr1Alt-gr2Alt; maxdiffgrAlt=max(diffgrAlt(:)); mindiffgrAlt=min(diffgrAlt(:)); mingrAlt=min(gr1Alt(:)); maxgrAlt=max(gr1Alt(:)); mediangrAlt=median(gr1Alt(:)); mediandiffgrAlt=median(diffgrAlt(:));
 
%SAG
%Max ve min Diff bul
rg1Sag=rg(BorderY:Border2Y,sizeX-Border-1:sizeX-1); rg2Sag=rg(BorderY:Border2Y,sizeX-Border:sizeX); diffrgSag=rg1Sag-rg2Sag; maxdiffrgSag=max(diffrgSag(:)); mindiffrgSag=min(diffrgSag(:)); minrgSag=min(rg1Sag(:)); maxrgSag=max(rg1Sag(:)); medianrgSag=median(rg1Sag(:)); mediandiffrgSag=median(diffrgSag(:));
rb1Sag=rb(BorderY:Border2Y,sizeX-Border-1:sizeX-1); rb2Sag=rb(BorderY:Border2Y,sizeX-Border:sizeX); diffrbSag=rb1Sag-rb2Sag; maxdiffrbSag=max(diffrbSag(:)); mindiffrbSag=min(diffrbSag(:)); minrbSag=min(rb1Sag(:)); maxrbSag=max(rb1Sag(:)); medianrbSag=median(rb1Sag(:)); mediandiffrbSag=median(diffrbSag(:)); 
 
br1Sag=br(BorderY:Border2Y,sizeX-Border-1:sizeX-1); br2Sag=br(BorderY:Border2Y,sizeX-Border:sizeX); diffbrSag=br1Sag-br2Sag; maxdiffbrSag=max(diffbrSag(:)); mindiffbrSag=min(diffbrSag(:)); minbrSag=min(br1Sag(:)); maxbrSag=max(br1Sag(:)); medianbrSag=median(br1Sag(:)); mediandiffbrSag=median(diffbrSag(:));
bg1Sag=bg(BorderY:Border2Y,sizeX-Border-1:sizeX-1); bg2Sag=bg(BorderY:Border2Y,sizeX-Border:sizeX); diffbgSag=bg1Sag-bg2Sag; maxdiffbgSag=max(diffbgSag(:)); mindiffbgSag=min(diffbgSag(:)); minbgSag=min(bg1Sag(:)); maxbgSag=max(bg1Sag(:)); medianbgSag=median(bg1Sag(:)); mediandiffbgSag=median(diffbgSag(:));
 
gb1Sag=gb(BorderY:Border2Y,sizeX-Border-1:sizeX-1); gb2Sag=gb(BorderY:Border2Y,sizeX-Border:sizeX); diffgbSag=gb1Sag-gb2Sag; maxdiffgbSag=max(diffgbSag(:)); mindiffgbSag=min(diffgbSag(:)); mingbSag=min(gb1Sag(:)); maxgbSag=max(gb1Sag(:)); mediangbSag=median(gb1Sag(:)); mediandiffgbSag=median(diffgbSag(:));
gr1Sag=gr(BorderY:Border2Y,sizeX-Border-1:sizeX-1); gr2Sag=gr(BorderY:Border2Y,sizeX-Border:sizeX); diffgrSag=gr1Sag-gr2Sag; maxdiffgrSag=max(diffgrSag(:)); mindiffgrSag=min(diffgrSag(:)); mingrSag=min(gr1Sag(:)); maxgrSag=max(gr1Sag(:)); mediangrSag=median(gr1Sag(:)); mediandiffgrSag=min(diffgrSag(:));
 
 
%?uan Hepsinin - rg, rb..- max ve min difleri directionlar ayri olarak var 
%birle?tirilmesi ve yeni diffmin ve diffmax bulunmas? laz?m
%min
mindiffrgarr=[mindiffrgUst mindiffrgSol mindiffrgAlt mindiffrgSag]; mindiffrg=min(mindiffrgarr(:));
mindiffrbarr=[mindiffrbUst mindiffrbSol mindiffrbAlt mindiffrbSag]; mindiffrb=min(mindiffrbarr(:));
mindiffbrarr=[mindiffbrUst mindiffbrSol mindiffbrAlt mindiffbrSag]; mindiffbr=min(mindiffbrarr(:));
mindiffbgarr=[mindiffbgUst mindiffbgSol mindiffbgAlt mindiffbgSag]; mindiffbg=min(mindiffbgarr(:));
mindiffgbarr=[mindiffgbUst mindiffgbSol mindiffgbAlt mindiffgbSag]; mindiffgb=min(mindiffgbarr(:));
mindiffgrarr=[mindiffgrUst mindiffgrSol mindiffgrAlt mindiffgrSag]; mindiffgr=min(mindiffgrarr(:));
%max
maxdiffrgarr=[maxdiffrgUst maxdiffrgSol maxdiffrgAlt maxdiffrgSag]; maxdiffrg=max(maxdiffrgarr(:));
maxdiffrbarr=[maxdiffrbUst maxdiffrbSol maxdiffrbAlt maxdiffrbSag]; maxdiffrb=max(maxdiffrbarr(:));
maxdiffbrarr=[maxdiffbrUst maxdiffbrSol maxdiffbrAlt maxdiffbrSag]; maxdiffbr=max(maxdiffbrarr(:));
maxdiffbgarr=[maxdiffbgUst maxdiffbgSol maxdiffbgAlt maxdiffbgSag]; maxdiffbg=max(maxdiffbgarr(:));
maxdiffgbarr=[maxdiffgbUst maxdiffgbSol maxdiffgbAlt maxdiffgbSag]; maxdiffgb=max(maxdiffgbarr(:));
maxdiffgrarr=[maxdiffgrUst maxdiffgrSol maxdiffgrAlt maxdiffgrSag]; maxdiffgr=max(maxdiffgrarr(:));
 
 
%?imdi al?nan parçalar?n max ve minlerini bulmak laz?m 
%birle?tirip tekrar bulaca??z.
 
 
minrgarr=[minrgUst minrgSol minrgAlt minrgSag]; minrg=min(minrgarr(:));
minrbarr=[minrbUst minrbSol minrbAlt minrbSag]; minrb=min(minrbarr(:));
minbrarr=[minbrUst minbrSol minbrAlt minbrSag]; minbr=min(minbrarr(:));
minbgarr=[minbgUst minbgSol minbgAlt minbgSag]; minbg=min(minbgarr(:));
mingbarr=[mingbUst mingbSol mingbAlt mingbSag]; mingb=min(mingbarr(:));
mingrarr=[mingrUst mingrSol mingrAlt mingrSag]; mingr=min(mingrarr(:));
%max
maxrgarr=[maxrgUst maxrgSol maxrgAlt maxrgSag]; maxrg=max(maxrgarr(:));
maxrbarr=[maxrbUst maxrbSol maxrbAlt maxrbSag]; maxrb=max(maxrbarr(:));
maxbrarr=[maxbrUst maxbrSol maxbrAlt maxbrSag]; maxbr=max(maxbrarr(:));
maxbgarr=[maxbgUst maxbgSol maxbgAlt maxbgSag]; maxbg=max(maxbgarr(:));
maxgbarr=[maxgbUst maxgbSol maxgbAlt maxgbSag]; maxgb=max(maxgbarr(:));
maxgrarr=[maxgrUst maxgrSol maxgrAlt maxgrSag]; maxgr=max(maxgrarr(:));
 
 
sh=0.0725;
diffrgarr=[maxdiffrg (-1)*mindiffrg]; diffrg=max(diffrgarr(:))+sh;
diffrbarr=[maxdiffrb (-1)*mindiffrb]; diffrb=max(diffrbarr(:))+sh;
diffbrarr=[maxdiffbr (-1)*mindiffbr]; diffbr=max(diffbrarr(:))+sh;
diffbgarr=[maxdiffbg (-1)*mindiffbg]; diffbg=max(diffbgarr(:))+sh;
diffgbarr=[maxdiffgb (-1)*mindiffgb]; diffgb=max(diffgbarr(:))+sh;
diffgrarr=[maxdiffgr (-1)*mindiffgr]; diffgr=max(diffgrarr(:))+sh;
 
%Thresholdlar belirlenir
threshrg1=diffrg+maxrg; threshrg2=((-1)*diffrg)+minrg;
threshrb1=diffrb+maxrb; threshrb2=((-1)*diffrb)+minrb;
threshbr1=diffbr+maxbr; threshbr2=((-1)*diffbr)+minbr;
threshbg1=diffbg+maxbg; threshbg2=((-1)*diffbg)+minbg;
threshgb1=diffgb+maxgb; threshgb2=((-1)*diffgb)+mingb;
threshgr1=diffgr+maxgr; threshgr2=((-1)*diffgr)+mingr;
 
 
 
 
TRG= rg>threshrg1 | rg<threshrg2; rg(TRG)=255;rg(~TRG)=0;
TRB= rb>threshrb1 | rb<threshrb2; rb(TRB)=255;rb(~TRB)=0;
 
TBR= br>threshbr1 | br<threshbr2; br(TBR)=255;br(~TBR)=0;
TBG= bg>threshbg1 | bg<threshbg2; bg(TBG)=255;bg(~TBG)=0;
 
TGB= gb>threshgb1 | gb<threshgb2; gb(TGB)=255;gb(~TGB)=0;
TGR= gr>threshgr1 | gr<threshgr2; gr(TGR)=255;gr(~TGR)=0;
 
final=(rg&rb)|(br&bg)|(gb&gr); 
 
Final(~final)=0;Final(final)=255; Final(:,sizeX)=0;
 
 
 
 
 
 
 
 
 
rg2=r./g;rb2=r./b;
br2=b./r;bg2=b./g;
gb2=g./b;gr2=g./r;
 
 
medRG=median(rg2(final)); medRB=median(rb2(final));
medBR=median(br2(final)); medBG=median(bg2(final));
medGB=median(gb2(final)); medGR=median(gr2(final));
 
maxRG=max(rg2(final)); maxRB=max(rb2(final));
maxBR=max(br2(final)); maxBG=max(bg2(final));
maxGB=max(gb2(final)); maxGR=max(gr2(final));
 
minRG=min(rg2(final)); minRB=min(rb2(final));
minBR=min(br2(final)); minBG=min(bg2(final));
minGB=min(gb2(final)); minGR=min(gr2(final));
 
diffRG1=maxRG-medRG; diffRB1=maxRB-medRB;
diffBR1=maxBR-medBR; diffBG1=maxBG-medBG;
diffGB1=maxGB-medGB; diffGR1=maxGR-medGR;
 
diffRG2=(-1)*minRG+medRG; diffRB2=(-1)*minRB+medRB;
diffBR2=(-1)*minBR+medBR; diffBG2=(-1)*minBG+medBG;
diffGB2=(-1)*minGB+medGB; diffGR2=(-1)*minGR+medGR;
 
x=find(Final==255);
y=find(Final==0);
if isempty(diffRG1) || isempty(diffRB1) || isempty(diffBR1) || isempty(diffBG1) || isempty(diffGB1) || isempty(diffGR1) || isempty(diffRG2) || isempty(diffRB2) || isempty(diffBR2) || isempty(diffBG2) || isempty(diffGB2) || isempty(diffGR2) || (size(x,1)<5000)|| (size(y,1)<20000) || size(x,1) > 4*size(y,1)
r = I(:,:,1);
g = I(:,:,2);
b = I(:,:,3);
%im to double
r=double(r);
g=double(g);
b=double(b);
 
%First segmentation
t= ((b./g)<0.95) & ((r./g)<0.95);
 
 
%Second segmentation
Igray=b;
 
T1=0;
T2=mean2(Igray);
T=(T1+T2)/2;
 
    done=false;
    while ~done
 
        %Background pixels
        r1=Igray<=T;
        %Object pixels
        r2=Igray>T;
 
        %New threshold
        T1=mean(Igray(r1));
        T2=mean(Igray(r2));
        Tnew=(T1+T2)/2;
 
        %Converge 
        done=abs(Tnew-T)<1;
        T=Tnew;
    end
  
 
%Combine pixels
r1=r1|t;
    
%If pixel is smaller than threshold
Igray(r1)=255;
%If pixel is bigger than threshold
Igray(~r1)=0;
 
%Morphological operations-fill
se = strel('square',5);  
Igray = imclose(Igray,se);
 
%Result
seg=Igray;
else
 
    
consrg1=diffRG1^2-diffRG1^4; consrg2=diffRG2^2-diffRG2^4;
consrb1=diffRB1^2-diffRB1^4; consrb2=diffRB2^2-diffRB2^4;
consbr1=diffBR1^2-diffBR1^4; consbr2=diffBR2^2-diffBR2^4;
consbg1=diffBG1^2-diffBG1^4; consbg2=diffBG2^2-diffBG2^4;
consgb1=diffGB1^2-diffGB1^4; consgb2=diffGB2^2-diffGB2^4;
consgr1=diffGR1^2-diffGR1^4; consgr2=diffGR2^2-diffGR2^4;
 
 
 
%if shadow
if diffRG1<diffRG2
    rg2=rg2<(medRG-diffRG1)+consrg1;
else
    rg2=rg2>(medRG+diffRG2)-consrg2;
end
%
if diffRB1<diffRB2
    rb2=rb2<(medRB-diffRB1)+consrb1;
else
    rb2=rb2>(medRB+diffRB2)-consrb2;
end
%
if diffBR1<diffBR2
    br2=br2<(medBR-diffBR1)+consbr1;
else
    br2=br2>(medBR+diffBR2)-consbr2;
end
%
if diffBG1<diffBG2
    bg2=bg2<(medBG-diffBG1)+consbg1;
else
    bg2=bg2>(medBG+diffBG2)-consbg2;
end
%
if diffGB1<diffGB2
    gb2=gb2<(medRG-diffGB1)+consgb1;
else
    gb2=gb2>(medGB+diffGB2)-consgb2;
end
%
if diffGR1<diffGR2
    gr2=gr2<(medGR-diffGR1)+consgr1;
else
    gr2=gr2>(medGR+diffGR2)-consgr2;
end
%
 
final2=(rg2&rb2)|(br2&bg2)|(gb2&gr2); 
 
Final(final)=255;   
Final(final2)=0;
Final(~final)=0;
 
 
 seg=Final;
end


e1 = edge(I(:,:,1), 'canny', 0.25, 2);
AX=seg==255;
denk=AX|e1;
seg(denk)=255;
seg(~denk)=0;

se = strel('disk',5); 

seg=imclose(seg,se);
seg=imfill(seg);
se = strel('disk',5); 

seg = imerode(seg,se);
se = strel('disk',2); 
seg=imdilate(seg,se);


sizeof0=size(find(seg==0),1);
sizeof255=size(find(seg==255),1);
if((sizeof0<20000||sizeof255<7500))

    %Image Size
    [nrows,ncols,nchans]=size(INEW);

    % Gaussian Filtering
    h = fspecial('gaussian', 5, 5);
    I_gauss = imfilter(INEW,h);

    % Influence of color data
    THETA = 1;
    % Influence of edge data
    BETA = 0.3;

    % Edge Data of RGB Image
    e1 = edge(INEW(:,:,1), 'canny', 0.25, 2);
   
    e2 = edge(INEW(:,:,2), 'canny', 0.25, 2);
    e3 = edge(INEW(:,:,3), 'canny', 0.25, 2);
    edgeData = double(reshape(e1+e2+e3, nrows*ncols,1));



    % Changing to lab color space
    cform = makecform('srgb2lab');
    lab = applycform(I_gauss,cform);


    % Calculating variance of L and z by reshaping matrix
    var_L = sqrt(var(double(reshape(lab(:,:,1), 1, ncols*nrows))));

    % Normalized L*a*b* and xyz information
    lab_ = (1/var_L)*double(lab);


    % Constructing feature vector
    feature_vec = THETA*double([reshape(lab_,nrows*ncols,3)...
                  BETA*edgeData]);


    opts = statset('MaxIter',300);
    % do k-means with 2 clusters
    idx= kmeans(feature_vec, 2,...
                                   'replicate', 3,...
                                   'EmptyAction', 'drop',...
                                   'Options',opts);


    % convert to black-and-white image where first cluster is chosen
    BWimg=reshape(idx==1,nrows,ncols);

    X=BWimg==1;

    J=rgb2gray(INEW);
    J(X)=255;
    J(~X)=0;
    se = strel('disk',5); 
    J = imclose(J,se);
    J = imerode(J,se);
    J = imdilate(J,se);

    x1=J(:,35);
    x2=J(:,size(J,2))-35;
    x3=J(35,:);
    x4=J(size(J,1),:)-35;


    medx1=median(x1);
    medx2=median(x2);
    medx3=median(x3);
    medx4=median(x4);

    if((medx1>210&&medx2>210&&medx3>210)||(medx2>210&&medx3>210&&medx4>210)||(medx1>210&&medx2>210&&medx4>210)||(medx1>210&&medx3>210&&medx4>210) )

        Y=J==255;
        J(Y)=0;
        J(~Y)=255;

    end
    
    
    
    seg=J;
 
end

se = strel('disk',7); 
seg = imerode(seg,se);
seg = imdilate(seg,se);

se = strel('disk',3); 

seg = imdilate(seg,se);
 
end