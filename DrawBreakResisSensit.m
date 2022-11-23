function [] = DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,ModelId,Fd,Resist_Sensit)
   
     FigH = Plot_US_R_BreakResSensit1(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,Fd,Resist_Sensit,ModelId);   
     pause(5);
    
end

