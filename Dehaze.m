%Start the program
%Let us initialize the AHR
close all;
clear;
i=imread('1234.bmp');

%resize the image.
i=imresize(i,[720 1280]);
I=i;

%Dark Channel Prior Haze Remove Algorithm 
H=fspecial('average',[8 1]);
source=I;
dark=min(source,[],3);

dark2=imfilter(dark,H,'replicate');
% dark2=dark;
sumdark=sum(dark(:));
maxdark=max(dark(:));
maxl=max(source(:));
meandark=uint8(uint32(double(sumdark)*1.3/1280/720));
Tempp=min(230,meandark);
L=uint8(bitshift(uint16(dark2)*uint16(Tempp),-8));
L=min(dark,L);
imshow(L)
A=maxl/2+maxdark/2;
%Render A temp hazeremove image;
zzzzz=uint8(double(source-L)./double(1.0-double(L)/double(A)));
%Show original image
%imshow(I);
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
I=rgb2gray(I);
%!!!!Some Coefficients
numBins = 256;
normClipLimit = 0.005;
numTiles = [10 10];
fullRange=[0 255];
selectedRange  = fullRange;
dimI = size(I);
dimTile = dimI ./ numTiles;
%CLAHE Algorithm
numPixInTile = prod(dimTile);
minClipLimit = ceil(numPixInTile/numBins);
clipLimit =  ceil(normClipLimit*numPixInTile+numPixInTile/256);
[atileMappings parameters] = makeTileMappings(R,G,B, numTiles, dimTile, numBins, clipLimit,selectedRange, fullRange,I);
%Synthesize the output image based on the individual tile mappings. 
out(:,:,1) = makeClaheImage(R,atileMappings, numTiles, selectedRange, numBins,dimTile);
out(:,:,2) = makeClaheImage(G, atileMappings, numTiles, selectedRange, numBins,dimTile);
out(:,:,3) = makeClaheImage(B, atileMappings, numTiles, selectedRange, numBins,dimTile);
%figure;
%show CLAHE Image
%imshow(out); 
%Enlighten the DP image
ligt=min(mean(out(:)),180);
uo=ligt/mean(zzzzz(:));
zzzzz=uint8(uo*double(source-L)./double(1.0-double(L)/double(A)));

%figure;
%show dp image
%imshow(zzzzz);
%combine image
out2=uint8((uint16(out)*parameters(1)+uint16(zzzzz)*parameters(2))/8);

figure;
% show the combined image
imshow(out2);
imwrite(out2,'1234.jpg')




function [atileMappings parameters]= makeTileMappings(R,G,B, numTiles, dimTile, numBins, clipLimit,  selectedRange, fullRange,I)
numPixInTile = prod(dimTile);
meanmI=mean(I(:));
meanI=mean(I(numTiles(1):numTiles(2)));
tileMappings = cell([numTiles 3]);
meanItile=zeros(numTiles(1),numTiles(2));
% extract and process each tile
imgCol = 1;
for col=1:numTiles(2),
  imgRow = 1;
  for row=1:numTiles(1),
    Rtile = R(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    Gtile = G(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    Btile = B(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    Itile=I(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    meanItile(row,col)=mean(Itile(:));
    RtileHist = imhist(Rtile, numBins); 
    GtileHist = imhist(Gtile, numBins); 
    BtileHist = imhist(Btile, numBins); 
    atileMapping = clipHistogram(RtileHist,GtileHist,BtileHist, clipLimit, numBins,Rtile,Gtile,Btile,numPixInTile);
    % assemble individual tile mappings by storing them in a cell array;
    atileMappings{row,col} = atileMapping;
    imgRow = imgRow + dimTile(1);    
  end
  imgCol = imgCol + dimTile(2); % move to the next column of tiles
end
vrIt =std(meanItile(:));
sprintf('vrit%d\n',vrIt);
vrIt=9.5*meanmI/80*meanmI/80+vrIt;
sprintf('meanmI%d\n',meanmI);
if vrIt>70
    vrIt=70;
    end
if vrIt<5;
    vrIt=5;
end
    tvrIt=(vrIt/70)*(vrIt/70);
    parameters(1)=tvrIt*4.0+2.0;
    parameters(2)=4.0-tvrIt*4.0+2.0;
end



% This function clips the histogram according to the clipLimit and
% redistributes clipped pixels across bins below the clipLimit
function imgHist = clipHistogram(rimgHist,gimgHist,bimgHist,clipLimit, numBins,R,G,B,numPixInTile)
% total number of pixels overflowing clip limit in each bin
% this loop should speed up the operation by putting multiple pixels
% into the "obvious" places first
tempsum1=cumsum(rimgHist);
tempsum2=cumsum(gimgHist);
tempsum3=cumsum(bimgHist);
imgHist=floor((rimgHist+gimgHist+bimgHist)/3);
% total number of pixels overflowing clip limit in each bin
totalExcess = sum(max(imgHist - clipLimit,0));  
% clip the histogram and redistribute the excess pixels in each bin
avgBinIncr = floor(totalExcess/numBins);
upperLimit = clipLimit - avgBinIncr; % bins larger than this will be
% set to clipLimit
% this loop should speed up the operation by putting multiple pixels
% into the "obvious" places first
for k=1:numBins
  if imgHist(k) > clipLimit
    imgHist(k) = clipLimit;
  else
    if imgHist(k) > upperLimit % high bin count
      totalExcess = totalExcess - (clipLimit - imgHist(k));
      imgHist(k) = clipLimit;
    else
      totalExcess = totalExcess - avgBinIncr;
      imgHist(k) = imgHist(k) + avgBinIncr;      
    end
  end
end
% this loops redistributes the remaining pixels, one pixel at a time
k = 1;
while (totalExcess ~= 0)
  %keep increasing the step as fewer and fewer pixels remain for
  %the redistribution (spread them evenly)
  stepSize = max(floor(numBins/totalExcess),1);
  for m=k:stepSize:numBins
    if imgHist(m) < clipLimit
      imgHist(m) = imgHist(m)+1;
      totalExcess = totalExcess - 1; %reduce excess
      if totalExcess == 0
        break;
      end
    end
  end
k = k+1; %prevent from always placing the pixels in bin #1
  if k > numBins % start over if numBins was reached
    k = 1;
  end
end
histSum = cumsum(imgHist);
valSpread  = 255;
alpha=0.4;
  vmax = 1 - exp(-alpha);
  val = (vmax*histSum/numPixInTile);
  val(val>=1) = 1-eps;
  temp = -1/alpha*log(1-val);
  mapping = min(0+temp*valSpread,255);
  imgHist=mapping;
end



function claheI = makeClaheImage(I, tileMappings, numTiles, selectedRange,...
                                 numBins, dimTile)
%initialize the output image to zeros (preserve the class of the input image)
claheI = I;
claheI(:) = 0;
%compute the LUT for looking up original image values in the tile mappings,
%which we created earlier
aLut=0:255;
imgTileRow=1;
for k=1:numTiles(1)+1
  if k == 1  %special case: top row
    imgTileNumRows = dimTile(1)/2; %always divisible by 2 because of padding
    mapTileRows = [1 1];
  else 
    if k == numTiles(1)+1 %special case: bottom row      
      imgTileNumRows = dimTile(1)/2;
      mapTileRows = [numTiles(1) numTiles(1)];
    else %default values
      imgTileNumRows = dimTile(1); 
      mapTileRows = [k-1, k]; %[upperRow lowerRow]
    end
  end
  % loop over columns of the tileMappings cell array
  imgTileCol=1;
  for l=1:numTiles(2)+1
    if l == 1 %special case: left column
      imgTileNumCols = dimTile(2)/2;
      mapTileCols = [1, 1];
    else
      if l == numTiles(2)+1 % special case: right column
        imgTileNumCols = dimTile(2)/2;
        mapTileCols = [numTiles(2), numTiles(2)];
      else %default values
        imgTileNumCols = dimTile(2);
        mapTileCols = [l-1, l]; % right left
      end
    end
    % Extract four tile mappings
    ulMapTile = tileMappings{mapTileRows(1), mapTileCols(1)};
    urMapTile = tileMappings{mapTileRows(1), mapTileCols(2)};
    blMapTile = tileMappings{mapTileRows(2), mapTileCols(1)};
    brMapTile = tileMappings{mapTileRows(2), mapTileCols(2)};
    % Calculate the new greylevel assignments of pixels 
    % within a submatrix of the image specified by imgTileIdx. This 
    % is done by a bilinear interpolation between four different mappings 
    % in order to eliminate boundary artifacts.
    normFactor = imgTileNumRows*imgTileNumCols; %normalization factor  
    imgTileIdx = {imgTileRow:imgTileRow+imgTileNumRows-1, ...
                 imgTileCol:imgTileCol+imgTileNumCols-1};
  %  imgPixVals =uint8(fuzz(I(imgTileIdx{1},imgTileIdx{2}), aLut));
    imgPixVals =I(imgTileIdx{1},imgTileIdx{2});
    % calculate the weights used for linear interpolation between the
    % four mappings
    rowW = repmat((0:imgTileNumRows-1)',1,imgTileNumCols);
    colW = repmat(0:imgTileNumCols-1,imgTileNumRows,1);
    rowRevW = repmat((imgTileNumRows:-1:1)',1,imgTileNumCols);
    colRevW = repmat(imgTileNumCols:-1:1,imgTileNumRows,1);
    
    claheI(imgTileIdx{1}, imgTileIdx{2}) = ...
        (rowRevW .* (colRevW .* double(findt(imgPixVals,ulMapTile)) + ...
                     colW    .* double(findt(imgPixVals,urMapTile)))+ ...
         rowW    .* (colRevW .* double(findt(imgPixVals,blMapTile)) + ...
                     colW    .* double(findt(imgPixVals,brMapTile))))...
        /normFactor;
    
    imgTileCol = imgTileCol + imgTileNumCols;    
  end %over tile cols
  imgTileRow = imgTileRow + imgTileNumRows;
end %over tile rows
claheI=im2uint8(claheI);
end



%just look up table for 256bins
function zzz=findt(sss,zz)
zzz=zz(uint8(sss)+1);
end