function   promedio= avgroi(I,mask)
% promedio= avgroi(I,mask)
[~,~,s]=size(I);
promedio=zeros(s,1);
mask = logical(mask);
%%

for j=1:1:s
    S=I(:,:,j);
    promedio(j)=mean(  removeNaN( S(mask)));
end
%%