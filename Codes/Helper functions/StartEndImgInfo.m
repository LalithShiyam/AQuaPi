function [StartSlice,EndSlice]= StartEndImgInfo(Volume)

for lp=1:size(Volume,3)
    if nnz(Volume(:,:,lp))==0
       continue
    else
        StartSlice=lp;
        break
    end
end

for lp=size(Volume,3):-1:1
    if nnz(Volume(:,:,lp))==0
        continue
    else
        EndSlice=lp;
        break
    end
end

end
