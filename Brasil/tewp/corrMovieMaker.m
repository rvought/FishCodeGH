

    writerObj = VideoWriter('tensoranimation-f22s-1200s.avi');
    writerObj.FrameRate = 15;
    open(writerObj);


for j=1:5:1000 
    asdf = bTensor(f, 22, [j, j+60]); 
    mframe= getframe(gcf);
    writeVideo(writerObj,mframe);    

    pause(0.01); 
end
    
    
close(writerObj);
