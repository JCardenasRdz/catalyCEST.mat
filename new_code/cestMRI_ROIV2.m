%% CestMRI_ROI

function   [CESTfit_ROI] = cestMRI_ROIV2(cest_mri_prescan,cest_mri_postscan,ppm,mask,Method)
%% 1) Apply gaussian filter to raw CEST MRI data
gaussFilter = fspecial('gaussian',[3 3],1);
cest_mri_prescan= imfilter(cest_mri_prescan,gaussFilter,'same');
cest_mri_postscan= imfilter(cest_mri_postscan,gaussFilter,'same');
%% 2) Define variables from method structure:

x0=Method.x0;       Npools=Method.Npools;   point1=Method.range(1);
                                            point2=Method.range(2);


%% 2) Get dimensions of acidoCEST data
[rows,cols,~]=size(cest_mri_postscan);


%% 3) Define options for NLSQ Method
% Fitting Method
options = optimoptions('lsqcurvefit');  
% Max. number of iterations
options.MaxIter=20E3; 
% Max. number of func. evaluations
options.MaxFunEvals=20E3;                
% Tolerance for NLSQ method
options.TolX=1e-10;             options.TolFun=1e-10;
% Turn iterative display off
options.Display='off';

%%  4) Reshape mask to vectorize and preallocate maps
% Transform
mask=reshape(mask,rows*cols,[]);
% preallocation
MAPS.rsq=nan(size(mask));  
MAPS.cest=nan(rows*cols,Npools);        % a 2D matrix     
MAPS.ppm=MAPS.cest;                     % a 2D matrix 


%%  5)  Find indexes for pixeles to be analyzed
% note: All pixels outside ROI are equalt to zero in the binary mask
indices= find(mask); 

%%  6)  preallocate CESTfit cell
CESTfit_ROI=cell(size(indices));

%%  7)  Reshape CEST MRI 3D matrix into a 2D matrix
cest_mri_prescan=reshape(cest_mri_prescan,rows*cols,[]);
cest_mri_prescan=cest_mri_prescan(:,point1:point2);
cest_mri_postscan=reshape(cest_mri_postscan,rows*cols,[]);
cest_mri_postscan=cest_mri_postscan(:,point1:point2);
%%  8)  Define paramters needed to allocate maps
offset_index=(Npools*3)-(Npools-1):(Npools*3);
ppm=ppm(point1:point2)';   

%% 9) Calculate Z spectrums and Difference and perform Lorentzian Line Fitting
matrix_prescan=[];
for n=1:length(indices);
    Signal_prescan=cest_mri_prescan(indices(n),:)';
    matrix_prescan=cat(2,matrix_prescan, Signal_prescan);
end
Zspectrum_prescan=mean(matrix_prescan,2);
Zspectrum_prescan=Zspectrum_prescan./Zspectrum_prescan(1);

matrix_postscan=[];
for n=1:length(indices);
    Signal_postscan=cest_mri_postscan(indices(n),:)';
    matrix_postscan=cat(2,matrix_postscan, Signal_postscan);
end
Zspectrum_postscan=mean(matrix_postscan,2);
Zspectrum_postscan=Zspectrum_postscan./Zspectrum_postscan(1);

[~,i]=min(Zspectrum_postscan); ppmadj=ppm-ppm(i);


%% modified by Julio on 08-Aug-2017
%Difference=Zspectrum_prescan-Zspectrum_postscan;

%% 9.4) Extrapolate and fit Zspectrum to a N-lorentzian model
%ppmlong=linspace(min(ppmadj),max(ppmadj),1000); ppmlong=ppmlong';

% cubic spline before difference
Difference= Zspectrum_postscan - Zspectrum_prescan;

%[beta,~,R]=lsqcurvefit(@lorentzian,x0(:,2),ppmlong,Zlong,x0(:,1),x0(:,3),options);
    
[beta,~,R]=lsqcurvefit(@lorentzian,x0(:,2),ppmadj,Difference,x0(:,1),x0(:,3),options);
   
%% 9.5) Allocate fitting results for voxel Q 
CESTfit_ROI=cell(1);
CESTfit_ROI{1}.cfit.pars=beta;
CESTfit_ROI{1}.cfit.residuals=R;

   
    [Lsum,L]=lorentzian(beta,ppmadj);
            CESTfit_ROI{1}.Lsum=Lsum;
            CESTfit_ROI{1}.cfitall=L;
            CESTfit_ROI{1}.Zspectrum_prescan=Zspectrum_prescan;
            CESTfit_ROI{1}.Zspectrum_postscan=Zspectrum_postscan;
            CESTfit_ROI{1}.Difference=Difference;
            CESTfit_ROI{1}.ppmadj=ppmadj';
%% 9.6) Calculate FWHM and CEST amplitudes
x = ppmadj;
y= CESTfit_ROI{1}.cfitall(:,1);

% Julio on 08-Aug-2017

CESTfit_ROI{1}.rsq = rsquare(CESTfit_ROI{1}.Difference,CESTfit_ROI{1}.Lsum);

z= CESTfit_ROI{1}.cfitall(:,2);
maxCEST1 = max(y); 
maxCEST2 = max(z);
f = find(y==maxCEST1); 
cp = x(f);
y1= y./maxCEST1;
ydatawr(:,1) = y1;
ydatawr(:,2) = x;
newFit1=find(x>= cp);
newFit2=find(x < cp);
ydatawr2 = ydatawr(min(newFit1):max(newFit1),:);
ydatawr3 = ydatawr(min(newFit2):max(newFit2),:);
sp1 = spline(ydatawr2(:,1),ydatawr2(:,2),0.5);
sp2 = spline(ydatawr3(:,1),ydatawr3(:,2),0.5);
Fullw = sp1-sp2;

CESTfit_ROI{1}.fullwidth=Fullw;
CESTfit_ROI{1}.CEST1amp=maxCEST1;
CESTfit_ROI{1}.CEST2amp=maxCEST2;
end