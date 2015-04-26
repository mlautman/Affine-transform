function [mask,P] = RandomWalkerSeg(datavector,Beta)
%Random Walker Algorithm for Segmentation Challenge
%Image is the image as a matrix of intensities
%Beta is a free parameter used in the weighting function

%% Define the Lattice Structure based on the size of the Image
image=datavector.img;
seed=datavector.seed;


% Obtain the dimensions of the image
[xdim,ydim,o]=size(image);

% Create a grid of those dimensions starting at 0
%xrange=0:(xdim-1);
%yrange=0:(ydim-1);
%[xgrid, ygrid]=meshgrid(xrange,yrange);
%gridpoints=[xgrid(:),ygrid(:)];
Npoints=xdim*ydim;

% Define the edges connecting the points
edges=edges4connected(xdim,ydim);
%edges=[[(1:Npoints);(1:Npoints)'],[((1:Npoints)+1)';(1:Npoints)'+xdim]];

% Remove the outside edges and edges that do not match the size of the image
%toberemoved=find((edges(:,1)>Npoints)|(edges(:,1)<1)|(edges(:,2)>Npoints)|(edges(:,2)<1));
%edges([toberemoved;(xdim:xdim:((ydim-1)*xdim))'],:)=[];

%% Define the Weights from the Random Walker formulation

%Calculate intensity differences (Ii-Ij)^2
%Idiff=sqrt(sum((Image(edges(:,1),:)-Image(edges(:,2),:)).^2,2));
Idiff=abs(sum((image(edges(:,1))-image(edges(:,2))),2));

%Find the maximum difference value to use as the normalization constant rho
rho=max(Idiff);

%Define a small constant epsilon to improve calculations
eps=10^-6; %Based on paper

%Define the weight function as described in the paper
wij=exp(Beta/rho*(Idiff))+eps;

%clear Beta 
clear rho
clear Idiff
clear eps
%% Calculate the sparse Laplacian matrix L with elements Lm, B, B', and Lu

% Calculate a sparse 2D matrix of the weights, which will place weights at
% specific points where there are adjacent nodes and ignore zeros
% (non-adjacent nodes based on the formulation for Lij) for optimization
% purposes
sparseweights=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[wij;wij],Npoints,Npoints);

clear edges 
clear Npoints
clear wij

% Calculate L with Lij=di (sum of the wij terms along i minus wij) for i=j and -wij for adjacent points
% The diag function in Matlab assigns values to i=j, so by assigning the
% diagonal values as the sum of the sparseweights, one can subtract the
% sparseweights to the entire diagonal matrix to get di=sum(wij,i)-wij and -wij at adjacent off-diagonals 

L=diag(sum(sparseweights));
L=L-sparseweights;

clear sparseweights

%% Assign labels from the seed

% Find the number of seeds, which equals the number of labels
sparseseed=sparse(seed);
nseeds=nnz(sparseseed);
%nlabels=nseeds;
labels=unique(sparseseed);
labels=unique(sparseseed(find(sparseseed>0)));
nlabels=length(labels);

% Assign the label indices
labelindex=linspace(1,nlabels,nlabels);

%Let's try this instead--I think I Need to have the indices as the seed
%locations (the soods you were trying last night)
label1loc=find(sparseseed==1);
label2loc=find(sparseseed==2);
label3loc=find(sparseseed==3);
%label4loc=find(sparseseed==4);

labelindex=[label1loc(length(label1loc)/2),label2loc(length(label2loc)/2),label3loc(length(label3loc)/2)];

% Find the values (in order) of the labels
%labels=nonzeros(sparseseed)';

% Define the matrix M which has elements msj=1 for s=j and zero elsewhere,
% where s and j are indices in the matrix. This matrix is used to separate
% B with the labels, so it should be an nxn diagonal matrix with n=nseeds
v=ones(1,nlabels);
M=diag(v);
M=sparse(M);

clear v
clear sparseseed

%% Calculate BT*M (The Right Hand Side of the Equation Lu*x=BT*M). 
%BT=-L(columnns outside the label index evaluated at each label)
coloutsidelabel=nlabels+1:length(L);
coloutsidelabel=sparse(coloutsidelabel);
BT=-L(coloutsidelabel,labelindex);
BT=sparse(BT);
RHS=BT*M;

clear BT
%% Use Matlab's implicit system of equations solver to find the potentials matrix X

X=L(coloutsidelabel,coloutsidelabel)\RHS;

%% Create a mask by selecting the label for which each point has the maximum potential

% Combine M and X so that in the final calculation the seeds are labeled
% correctly (potential=probability=1 for each seed)
P=[M;X];

% Find the label to which each point should belong by taking the maximum
% along each row (dimension 2) and finding the index for that point. 
[maxvalue, assignedlabel]=max(P,[],2);

%Combine all of these labels into a max in the dimensions of the original
%image
mask=reshape(assignedlabel, [xdim,ydim]);
P=full(P);
P=reshape(P, [xdim,ydim,nlabels]);

end



