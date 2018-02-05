function [ areaWithParticles ] = initalizeParticles( numbPart, area )
%initalizeParticles initalizes particles evenly spaced
%   Initalize particles in the center of the sceen in x and y dimentations
%   provided by the number of xand y particles defined in inputs

%initalize constants
PartSeperation = 1e-9; %1 nm inital spacing
numy = 0; %inital posiotion for x
numx = 0; %inital posiotion for y

% %loop through the number x and rthen y particles, placing each in an inital position
% for indexX = 0:numbPartX - 1
%     numX = numX + 1;
%     x(numX) = x0 * PartSeperation - indexX; %* cos(PartAng);
%     
%     %x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
%     %AtomType(numbPartX) = Type;
% end
% 
% for indexY = 0:numbPartY - 1
%     numY = numY + 1;
%     y(numY) = y0 * PartSeperation - indexY;
% 
%     %y(numbPartX) = y0 * AtomSpacing - Seper * index * AtomSpacing * sin(PartAng);
%     %AtomType(numbPartX) = Type;
% end


%genreate numbPart trandome x locations, then numbPart randome y. 

 for index = 0:numbPart - 1
     PositionX(index) = rand; %generate randome number for x cordinate
     PositionY(index) = rand; %generate randome number for y cordinate
 end    
%      %verify particles are well seperated
%      if area(PositionX+PartSeperation)==1;
%          PositionX = rand; %re-generate x (too close right)
%      elseif  area(PositionY+PartSeperation)==1;
%          PositionY = rand; %re-generate y (too close right)
%      elseif area(PositionX-PartSeperation)==1;
%          PositionX = rand; %re-generate x (too close left)
%      elseif area(PositionY-PartSeperation)==1;
%          PositionY = rand; %re-generate y (too close left)
%      else         
%      area (PositionX, PositionY)= 1; %position particle in area (particle = 1, no particle =0)   
     
 %esure none are within partile seperaion requirement, else generate new
 %value replace that one,, recheck

for choose = 0:size(PositionX)-1
    for compare = 1:size(PositionX)-1

        %if location duplicated
        if PositionX(choose) == PositionX(compare) & PositionY(choose) == PositionY(compare)
             PositionX(choose) = rand; %re-generate x and y
             PositionY(choose) = rand;
             %dont increase compare index yet, we now need to recompare to each element
             compare=compare-1;        
             
        %if locations are too close right
        elseif PositionX(choose) ==  PositionX(compare)+ PartSeperation & PositionY(choose) == PositionY(compare)+ PartSeperation
             PositionX(choose) =  PositionX(choose) - PartSeperation;
             %dont increase compare index yet, we now need to recompare to each element
             compare=compare-1;
        %if locations are too close left
        elseif PositionX(choose) ==  PositionX(compare)- PartSeperation & PositionY(choose) == PositionY(compare)- PartSeperation
            PositionX(choose) =  PositionX(choose) + PartSeperation;
            %dont increase compare index yet, we now need to recompare to each element
            compare=compare-1;
            
       elseif PositionX(choose) ==  PositionX(compare)+ PartSeperation & PositionY(choose) == PositionY(compare)- PartSeperation
             PositionX(choose) =  PositionX(choose) - PartSeperation;
             %dont increase compare index yet, we now need to recompare to each element
             compare=compare-1;
        elseif PositionX(choose) ==  PositionX(compare)- PartSeperation & PositionY(choose) == PositionY(compare)+ PartSeperation
             PositionX(choose) =  PositionX(choose) - PartSeperation;
             %dont increase compare index yet, we now need to recompare to each element
             compare=compare-1;
    end
end

     
 end
 


