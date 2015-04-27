function [mask,P] = RandomWalkerSeg7(datavector)
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

%% Define a vector of Betas (range determined by testing) for beta selection
Betas=[90,95,100,111];
%Preallocate the initial dicefactor as 0. This needs to be outside of the
%loops to avoid extra calculations
diceface0=0;
for c=1:length(Betas)
    Beta=Betas(c);
    %% Define the Weights from the Random Walker formulation
    
    %Calculate intensity differences
    Idiff=abs(sum((image(edges(:,1))-image(edges(:,2))),2));
    
    %Find the maximum difference value to use as the normalization constant rho
    rho=max(Idiff)^2;
    
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
    L=sparse(diag(sum(W))-W);
    
    
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
    
    % Define the matrix M which has elements msj=1 for s=j and zero elsewhere,
    % where s and j are indices in the matrix. This matrix is used to separate
    % B with the labels, so it should be an nxn diagonal matrix with n=nseeds
    v=ones(1,nlabels);
    M=diag(v);
    M=sparse(M);
    
    clear v
    
    if nlabels==2
        for i=1:length(label1loc)
            for j=length(label2loc)
                seedindex=[label1loc(i),label2loc(j)];
                %seedindex=[label1loc(50),label2loc(50),label3loc(100)];
                
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
                
                X=sparse(L(coloutsidelabel,coloutsidelabel)\RHS);
                
                
                clear RHS
                
                %% Create a mask by selecting the label for which each point has the maximum potential
                
                % Combine M and X so that in the final calculation the seeds are labeled
                % correctly (potential=probability=1 for each seed)
                Pnew=[M;X];
                
                % Find the label to which each point should belong by taking the maximum
                % along each row (dimension 2) and finding the index for that point.
                [~, assignedlabel]=max(Pnew,[],2);
                
                %Combine all of these labels into a mask in the dimensions of the original
                %image
                masknew=reshape(assignedlabel, [xdim,ydim]);
                
                %When looking at the tumor file, we want to not look at the "other tissues"
                %area, so set those values to 0
                masknew(masknew>2)=0;
                gtr = gen_img_from_labels(datavector, datavector.seg);
                seg = gen_img_from_labels(datavector, masknew);
                [~,~,dice]=seg_eval(gtr,seg);
                
                %Calculate a determining factor made up of the average of
                %the dice scores. If the average dice score increases, it
                %keeps the mask from the iteration, otherwise it is not
                %saved
                dicefac=(dice(1,2)+dice(2,2));
                
                if dicefac>diceface0
                    mask=masknew;
                    P=Pnew;
                    %Update average dice0
                    diceface0=dicefac;
                end
            end
        end
        
    end
    
    if nlabels==3
        %seedindex=[label1loc(round(length(label1loc)/2)),label2loc(round(length(label2loc)/2)),label3loc(round(length(label3loc)/2))];
        for i=1:length(label1loc)
            for j=length(label2loc)
                for k=length(label3loc)
                    seedindex=[label1loc(i),label2loc(j),label3loc(k)];
                    %seedindex=[label1loc(50),label2loc(50),label3loc(100)];
                    
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
                    
                    X=sparse(L(coloutsidelabel,coloutsidelabel)\RHS);
                    
                    
                    clear RHS
                    
                    %% Create a mask by selecting the label for which each point has the maximum potential
                    
                    % Combine M and X so that in the final calculation the seeds are labeled
                    % correctly (potential=probability=1 for each seed)
                    Pnew=[M;X];
                    
                    % Find the label to which each point should belong by taking the maximum
                    % along each row (dimension 2) and finding the index for that point.
                    [~, assignedlabel]=max(Pnew,[],2);
                    
                    %Combine all of these labels into a mask in the dimensions of the original
                    %image
                    masknew=reshape(assignedlabel, [xdim,ydim]);
                    
                    %When looking at the tumor file, we want to not look at the "other tissues"
                    %area, so set those values to 0
                    masknew(masknew>2)=0;
                    gtr = gen_img_from_labels(datavector, datavector.seg);
                    seg = gen_img_from_labels(datavector, masknew);
                    [~,~,dice]=seg_eval(gtr,seg);
                    
                    %Calculate a determining factor made up of the average of
                    %the dice scores. If the average dice score increases, it
                    %keeps the mask from the iteration, otherwise it is not
                    %saved
                    dicefac=(dice(1,2)+dice(2,2));
                    
                    if dicefac>diceface0
                        mask=masknew;
                        P=Pnew;
                        %Update average dice0
                        diceface0=dicefac;
                    end
                    
                    
                end
            end
        end
    end
end
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
    title('Seeds')
    subplot(2,2,2)
    imagesc(image(:,:,1));colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    title('Image')
    subplot(2,2,3)
    imagesc(mask);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    title('Algorithm Segmentation')
    subplot(2,2,4)
    imagesc(datavector.seg);colormap('gray');hold on;plot(b1,a1,'g.','MarkerSize',24);plot(b2,a2,'b.','MarkerSize',24);plot(b3,a3,'r.','MarkerSize',24);
    title('True Segmentation')
end

end
