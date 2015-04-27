function [mask,P] = RandomWalkerSeg(datavector,Beta)
%Random Walker Algorithm for Segmentation Challenge
%Image is the image as a matrix of intensities
%Beta is a free parameter used in the weighting function
warning('off','all')
%% Define the Lattice Structure based on the size of the Image
image=datavector.img;

seed=datavector.seed;

% Obtain the dimensions of the image
[xdim,ydim,zdim]=size(image);
Npoints=xdim*ydim;

% Define the edges connecting the points
edges=edges4connected(xdim,ydim);


%% Define the Weights from the Random Walker formulation

%Calculate intensity differences
Idiff=abs(sum((image(edges(:,1))-image(edges(:,2))),2));

%Find the maximum difference value to use as the normalization constant rho
rho=max(Idiff);

%Define a small constant epsilon to improve calculations
eps=10^-6; %Based on paper

%Define the weight function as described in the paper
wij=exp(-Beta/rho*(Idiff))+eps;

clear Beta
clear rho
clear Idiff
clear eps
%% Calculate the sparse Laplacian matrix L with elements Lm, B, B', and Lu

%Create the Entire W matrix by shaping it along the edges and use to
%calculate L
W=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[wij;wij],Npoints,Npoints);
L=diag(sum(W))-W;


%% Assign labels from the seed

% Find the number of seeds, which equals the number of labels
sparseseed=sparse(seed);
labels=unique(sparseseed(find(sparseseed>0)));
labels=labels(labels>0);
nlabels=length(labels);

%Let's try this instead--I think I Need to have the indices as the seed
%locations (the soods you were trying last night)

label1loc=find(sparseseed==1);
label2loc=find(sparseseed==2);
label3loc=find(sparseseed==3);

clear sparseseed

if nlabels==2
    seedindex=[label1loc(length(label1loc)/2),label2loc(length(label2loc)/2)];
end

if nlabels==3
    seedindex=[label1loc(length(label1loc)/2),label2loc(length(label2loc)/2),label3loc(length(label3loc)/2)];
end
seedindex=[label1loc(50),label2loc(50),label3loc(100)];
% Define the matrix M which has elements msj=1 for s=j and zero elsewhere,
% where s and j are indices in the matrix. This matrix is used to separate
% B with the labels, so it should be an nxn diagonal matrix with n=nseeds
v=ones(1,nlabels);
M=diag(v);
M=sparse(M);

clear v
%% Calculate BT*M (The Right Hand Side of the Equation Lu*x=BT*M).
%BT=-L(columnns outside the label index evaluated at each label)

coloutsidelabel=1:length(L);
coloutsidelabel(seedindex)=[];
coloutsidelabel=sparse(coloutsidelabel);

BT=-L(coloutsidelabel,seedindex);
BT=sparse(BT);

%Calculate RHS
RHS=BT*M;

clear BT
%% Use Matlab's implicit system of equations solver to find the potentials matrix X

X=L(coloutsidelabel,coloutsidelabel)\RHS;

clear L
clear RHS

%% Create a mask by selecting the label for which each point has the maximum potential

% Combine M and X so that in the final calculation the seeds are labeled
% correctly (potential=probability=1 for each seed)
P=[M;X];

clear M

% Find the label to which each point should belong by taking the maximum
% along each row (dimension 2) and finding the index for that point.
[maxvalue, assignedlabel]=max(P,[],2);

%Combine all of these labels into a mask in the dimensions of the original
%image
mask=reshape(assignedlabel, [xdim,ydim]);

%When looking at the tumor file, we want to not look at the "other tissues"
%area, so set those values to 0
mask(mask>2)=0;

%Convert P from sparse to full
P=full(P);
P=reshape(P, [xdim,ydim,nlabels]);

clear figure
%%IMAGES
if nlabels==2
    [a1,b1]=ind2sub([240 240],seedindex(1));
    [a2,b2]=ind2sub([240 240],seedindex(2));
    subplot(2,2,1)
    imagesc(datavector.seed);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);
    subplot(2,2,2)
    imagesc(image(:,:,1));colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);
    subplot(2,2,3)
    imagesc(mask);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);
    subplot(2,2,4)
    imagesc(datavector.seg);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);
end
if nlabels==3
    [a1,b1]=ind2sub([240 240],seedindex(1));
    [a2,b2]=ind2sub([240 240],seedindex(2));
    [a3,b3]=ind2sub([240 240],seedindex(3));
    subplot(2,2,1)
    imagesc(datavector.seed);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    subplot(2,2,2)
    imagesc(image(:,:,1));colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    subplot(2,2,3)
    imagesc(mask);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    subplot(2,2,4)
    imagesc(datavector.seg);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
end

end
