%% CestMRI
% Estimate the CEST effect of an N number of pools by non-linear LSQ.
% This function assumes that a Z-spectrum can be modeled as the sumation of
% Lorentzian fucntions in the frequency domain.
%
%%  SYNTAX
% [CESTfit,MAPS,indices]=cestMRI(cest_mri,ppm,mask,Method)
% 
%%  INPUTS
%
% cest_mri: 3D Matrix.
%            where x=rows of image
%                  y=columns of images
%                  z=frequency offset of the transmitter
% ppm     : row vector of frequency offsets in ppm
% mask    : binary mask that defines the ROI to be analzyed
% 
% Method: Structure with the following fields.
%            Method.Npools= number of lorentzians functions;
%            Method.range = row vector with first and data points to be
%                           analyzed
%            Method.x0    =  3*Npools by 3*Npools matrix that contains
%                          Method.x0(:,1)=Lower Limits for NLSQ fitting
%                          Method.x0(:,2)=Initial guess for NLSQ fitting
%                          Method.x0(:,3)=Upper Limits for NLSQ fitting
%            each column of is composed of:
%            rows 1     to N:   Lorentzian amplitudes for pools one to N
%            rows N+1   to N*2: Lorentzian width for pools one to N
%            rows 2*N+1 to N*3: Lorentzian offsets for pools one to N
%%  OUTPUTS
%
% CESTfit{n}   Cell array of length n equal to the number of voxels
% analyzed. Each element contains the followint structure.

%   CESTfit{n}.cfit.pars:       Parameters obtained after NLSQ fitting
%   CESTfit{n}.cfit.residuals:  Residual of the fitting
%
%   CESTfit{n}.cfitsum:     Predicted Lorentzian.
%             .cfitall:     Matrix of predicted invidual Lorentzians
%             .Zspectrum:   Experimental Zspectrum (full length)
%             .ppmadj:      Adjusted offsets to center water
%
% MAPS  
% 
%      MAPS.rsq:   R square of the curve fitting
%           cest:  Matrix of cest amplitudes of size [rows,columns,Npools]
%           ppm:   Matrix of offsets  of size [rows,columns,Npools]
     
% indices: Vector indices for each pixel analyzed withint the ROI.
%
%%  EXAMPLE
%
%[CESTfit,MAPS,indices]=cestMRI(cest_mri,ppms,mask,Method);
% 
%% References:  
%   1)  Sheth VR, Li Y, Chen LQ, Howison CM, Flask CA, Pagel MD.
%       Measuring in vivo tumor pHe with CEST-FISP MRI. 
%       Magn. Reson. Med., 2012, 67:760?768.  PMID 22028287.
%
%   2)  Liu G, Li Y, Sheth VR, Pagel MD.
%       Imaging in vivo extracellular pH with a Single PARACEST MRI Contrast Agent.
%       Molecular Imaging, 2012, 11(1):47-57.  PMID 21651182.
%
%
%% Author
%   Julio Cárdenas-Rodríguez
%       The University of Arizona
%       cardenaj@email.arizona.edu
%       ver 1.0 March, 2014.

function   [CESTfit,MAPS,indices]=cestMRI_voxel(cest_mri_prescan,cest_mri_postscan,ppm,mask,Method)
%% 1) Apply gaussian filter to raw DCE MRI data
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
CESTfit=cell(size(indices));

%%  7)  Reshape CEST MRI 3D matrix into a 2D matrix
cest_mri_prescan=reshape(cest_mri_prescan,rows*cols,[]);
cest_mri_prescan=cest_mri_prescan(:,point1:point2);
cest_mri_postscan=reshape(cest_mri_postscan,rows*cols,[]);
cest_mri_postscan=cest_mri_postscan(:,point1:point2);
%%  8)  Define paramters needed to allocate maps
offset_index=(Npools*3)-(Npools-1):(Npools*3);
ppm=ppm(point1:point2)';   

%% 9) Perform Lorentzian Line Fitting only for Q voxels

for q=1:length(indices);
    %% 9.1) Get Signal for voxel Q and tranpose
    Signal_prescan= cest_mri_prescan(indices(q),:)';
    Signal_postscan= cest_mri_postscan(indices(q),:)';
    %% 9.2)Calculate Z spectrum
    %Signal=Signal-min(Signal);
    Zspectrum_prescan=Signal_prescan./Signal_prescan(1);
    Zspectrum_postscan=Signal_postscan./Signal_postscan(1);
   
    %% 9.3) Find Maximum and update ppm  to center Zspectrum
    [~,i]=min(Zspectrum_postscan); ppmadj=ppm-ppm(i);
    
    %% 9.4) Generate Difference Spectrum
    Difference=Zspectrum_prescan-Zspectrum_postscan;
    
    %% 9.4) Extrapolate and fit Zspectrum to a N-lorentzian model
ppmlong=linspace(min(ppmadj),max(ppmadj),1000); ppmlong=ppmlong';
Zlong = csaps(ppmadj,Difference,1,ppmlong);
%%
[beta,~,R,~,~,~,jacobian]=lsqcurvefit(@lorentzian,x0(:,2),ppmlong,Zlong,x0(:,1),x0(:,3),options);


%% modified by Julio on  03-Sep-2017 to include the confidence interval
ci = nlparci(beta,R,'jacobian',jacobian);

    %% 9.5) Allocate fitting results for voxel Q 
CESTfit{q}.cfit.pars=beta;
CESTfit{q}.cfit.residuals=R;
CESTfit{q}.cfit.ci=ci;

    [Lsum,L]=lorentzian(beta,ppmadj);
            CESTfit{q}.Lsum=Lsum;
            CESTfit{q}.cfitall=L;
            CESTfit{q}.Zspectrum_prescan=Zspectrum_prescan;
            CESTfit{q}.Zspectrum_postscan=Zspectrum_postscan;
            CESTfit{q}.Difference=Difference;
            CESTfit{q}.ppmadj=ppmadj';
            CESTfit{q}.S_prescan=Signal_prescan;
            CESTfit{q}.S_postscan=Signal_postscan;

rsq=rsquare(Difference,Lsum);
CESTfit{q}.cfit.rsq=rsq; 

   %% 9.6) Calculate FWHM and CEST amplitudes
x = ppmadj;
y= CESTfit{q}.cfitall(:,1);
z= CESTfit{q}.cfitall(:,2);
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

CESTfit{q}.fullwidth=Fullw;
CESTfit{q}.CEST1amp=maxCEST1;
CESTfit{q}.CEST2amp=maxCEST2;
    %%  9.7) Allocate parametric maps  
    MAPS.rsq(indices(q),1)=rsq;

    %% 9.8)
    % Allocate offsets and CEST effects
        for r=1:length(offset_index)
        MAPS.ppm(indices(q),r)=beta(offset_index(r));  
        end
    
        % Allocate CEST effects
        for r=1:length(offset_index)
        offset=beta(offset_index(r));  
        [~,l2]=lorentzian(beta,offset); 
        MAPS.cest(indices(q),r)=l2(r);
    
        end
    
end

%% 10) Reshapes MAPS back to square matrices

MAPS.cest=reshape(MAPS.cest,rows,cols,Npools); 
MAPS.rsq=reshape(MAPS.rsq,rows,cols);
MAPS.ppm=reshape(MAPS.ppm,rows,cols,Npools);

end




