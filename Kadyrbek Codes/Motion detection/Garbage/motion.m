function motion
v=videoinput('winvideo');
for r=1:30
    
    i=getsnapshot(v);
    i=rgb2gray(i);
    pause(0.2);
    j=getsnapshot(v);
    j=rgb2gray(j);
    z=imabsdiff(i,j);
    z=im2bw(z);
    [a b c]=size(z);
    x=0;
    for m=1:a
        for n=1:b
            if z(m,n)==1
                x=x+1;
            end
        end
    end    
    if(x>1000)
        disp('Motion Detected');
    else
        disp('No motion');
        
    end   
end
