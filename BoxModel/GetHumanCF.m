function [CF]=GetHumanCF(x,L,BM_flag,x0)

%The tonopic map for man adopted from Greenwood DD, A cochlear 
%frequency-position function for several species-29 years later

if BM_flag==1
   xo=1-x./L;
   A=165.4;
   a=2.1;
   k=0.88;
   CF=A*(10.^(a*xo)-k);
else
    xo=1-(x+x0)./L;
    A=165.4;
    a=2.1;
    k=0.88;
    CF=max(A*(10.^(a*xo)-k),5);
end

