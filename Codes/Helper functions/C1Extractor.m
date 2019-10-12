     %function [ C1SegmentICA ] = C1Extractor( ICA,FV,C1SliceEnd )
function [ C1SegmentICA ] = C1Extractor( FV,C1SliceEnd, OriginalVolume)

% Segments the C1 segment of the internal carotid artery
% 
% for lp=1:size(ICA,3)
%     [~,numberOfObject]=bwlabel(ICA(:,:,lp));
%     if numberOfObject > 1
%        C1Start=lp+1;
%        break
%     else
%         if numberOfObject==1 && lp==round(0.50*size(ICA,3))
%             C1Start=1;
%             break
%         else
%         end
%        continue
%     end
% end
% replace ica by fv


% patch added on 20-10-2016
tempFV=RemOverlay2(FV); %changed from removoerlay to removerlay2 : upgraded for vienna patients.
for lp=1:size(FV,3)
    se = strel('disk',3);
    temp(:,:,lp)=imerode(tempFV(:,:,lp),se); % it was ten before.
    [~,numberOfObject]=bwlabel(temp(:,:,lp));
    if numberOfObject>1
        BP=lp;
        break
    else
        if numberOfObject==1 && lp==round(0.5*size(FV,3))
            BP=1;
            break
        else
        end
    end  
end

C1SegmentICA=zeros(size(FV));
C1SegmentICA(:,:,BP+3:C1SliceEnd)=FV(:,:,BP+3:C1SliceEnd);

end

