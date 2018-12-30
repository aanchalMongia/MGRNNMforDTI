function [Interaction,S1,S2,Did,Tid]=getdata(dataset,similarity_type)
%GETDATA Summary of this function goes here

    load ('data/DtiData.mat', [dataset 'AdmatDGC_inl'] ,[dataset 'SimmatDC_inl'] , [dataset 'SimmatDG_inl'] )
    Interaction = eval([dataset 'AdmatDGC_inl']) ;
    S1=eval([dataset 'SimmatDC_inl']);
    S2= eval([dataset 'SimmatDG_inl'] );
    X=Interaction';
  
    
         for i=1:length(similarity_type)
          S1=S1+(1-pdist2(X,X,similarity_type{i}));
          S2=S2+(1-pdist2(X',X',similarity_type{i}));
         end

Did=[]; Tid=[];
end

