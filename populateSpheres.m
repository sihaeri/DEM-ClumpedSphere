function [ assembly, masses, totVolSph, mI ] = populateSpheres( fileName, nNodes, scaleFactor, smoothFact, writeTec, ...
                                                                writeLammpsTemp, lammpsType, writeStlScaled, rho, findSTLmI )
%=========================================================================%
%Copyright 2017 Sina Haeri                                                %
%This program is distributed under the terms of the GNU General Public    %  
%License                                                                  %
%=========================================================================%

%This function populates an STL volume with a number of spheres
%Spheres are non-overlapping and with pre-defined diameters. 
%This can on its own be used to generate a multi-sphere particle for DEM
%simulations but the number of spheres is often too large. Therefore, this will
%be used as an initial approximation to an optimum sphere filling algorithm
%and is not intended to be used directly. The next step is to use this
%initial assembly to generate an optimum cluster of larger spheres to
%represent the object for DEM simulations, the technique used is an improved 
%Greedy Algorithm. More details can found in Haeri (2017) Powder Technology
%94-104.
%--------------------------------------------------------------------------
%INPUTS:
%fileName: name of the STL file to be loaded. 
%nNodes: number of nodes to be placed along each axis (e.g. start with 50).  
%scaleFactor: scaling factor to normalise all coordinated in the STL file.
%scaleFactor: how much smoothing to apply to the surfce, this factor will
%be multiplied by the cell radius to smoothout the surface, if 1. is chosen
%then no smoothing and large number of spheres will be included in the
%final cluster to represent surface roughness (3-5 seems to be reasonable,
%see Haeri 2017 PowTech for more details. 
%writeTec: if = 1 then will write tecplot/paraview files of all the cells, 
%surface and internal cells and final cluster.
%writeLammpsTemp: Write a molecule template for lammps
%LammpsTypes: an integer value (the typew to be used in the lammps
%template, see lammps documentation for more details)
%writeStlScaled: if = 1, the scaled stl will be written
%rho: assumed density of the particle
%findSTLmI: 0 or 1, if 0 the Moment of Inertia for the STL is not
%calculated, this is very time consuming so default to 0, also only useful
%for comparison.
%--------------------------------------------------------------------------
%OUTPUTS:
%assembly: coordinate and the radius of the spheres in the final cluster
%masses: mass of the each sphere in the final cluster, this is calculated
%based on the number of cells exclusively inside each sphere in the cluster
% totVolSph: Volume of the cluster calculated based on the volume of the
% cells, this is not the volume of the original stl, that volume will be
% printed on the scree. 
% mI: moment of intertia of the cluster
%--------------------------------------------------------------------------

[vertices,faces,~,~] = stlRead(fileName);

vertices = vertices*scaleFactor;

%find a bounding box
xlow=min(vertices(:,1));
ylow=min(vertices(:,2));
zlow=min(vertices(:,3));

xhigh=max(vertices(:,1));
yhigh=max(vertices(:,2));
zhigh=max(vertices(:,3));

Lx = xhigh-xlow;
Ly = yhigh-ylow;
Lz = zhigh-zlow;

disp(['Bounding Box Dimensions Lx=', num2str(Lx), '; Ly=', num2str(Ly), ...
                                                  '; Lz=', num2str(Lz)]);

%Calculate some parameters
rSph = Lx/2/nNodes; %Radius of potential spherical cells

nRayY = floor(Ly/Lx*nNodes)+1; %Nomber of place holder lines in y-direction
nRayZ = floor(Lz/Lx*nNodes)+1; %Nomber of place holder lines in z-direction

%Determine the direction and initial origin of eah ray, direction is 
%arbitrarily set to +x, origins are distributed unifromely on x=xlow plane
%This should not be changed
initCoords = zeros(3,nRayY,nRayZ);
initDirs = zeros(3,nRayY,nRayZ);
for k = 1:nRayZ
    for j = 1:nRayY
        initCoords(:,j,k)=[xlow-1e13*eps(xlow); ylow+(j-1)*2*rSph; zlow+(k-1)*2*rSph]; 
        initDirs(:,j,k)=[Lx;0;0]; 
    end
end

%generate spheres along each ray (+x)
initSphCoords = zeros(nNodes+1,nRayY,nRayZ);
for k = 1:nRayZ
    for j = 1:nRayY
        for n = 1:nNodes+1
            initSphCoords(n,j,k)=xlow+(n-1)*2*rSph; 
        end
    end
end

%Just for testing a very simple crude Monte Carlo algorithm to calculte mI
%of the original STL, set findSTLmI to zero since this may take very long
%specially if the stl is not convex
[stlVol,~] = stlVolume(vertices',faces');
%stlMass = stlVol*rho;
if(findSTLmI) 
    [STLmI, STLCoM] = calculateSTLmI(vertices,faces,stlVol,rho,10000);
end

%For each ray determine the intersection
sz = size(faces);
nface= sz(1);

%intersections of each ray witht the collection of trinagles in +x
maxIntersectPerRay = 10; %Increase of crashed
allIntersects = zeros(maxIntersectPerRay, nRayY, nRayZ);
nIntersectRay = zeros(nRayY,nRayZ);

for k = 1:nRayZ
    for j = 1:nRayY
        direction = initDirs(:,j,k);
        origin    = initCoords(:,j,k);
        ncIntersect = 0;
        for facet=1:nface
            vert = vertices(faces(facet,:),:);
            [flag, ~, ~, t] = rayTriangleIntersection(origin, direction, vert(1,:)', vert(2,:)', vert(3,:)', eps(Lx));    
            if(flag)
                intersection = origin + t*direction;
                ncIntersect = ncIntersect+1;
                allIntersects(ncIntersect,j,k) = intersection(1);
            end
        end
        nIntersectRay(j,k) = ncIntersect;
        allIntersects(1:ncIntersect,j,k) = sort(allIntersects(1:ncIntersect,j,k));
    end
end

%bFlagSphere = 0: surface node, <7: inernal node, =7 
bFlagSphere = zeros(nNodes+1,nRayY,nRayZ);

for k = 1:nRayZ
    for j = 1:nRayY
        if(nIntersectRay(j,k)) %otherwise bFlagSphere is already set to 0
            for nn = 1:2:nIntersectRay(j,k)
                for n = 1:nNodes+1
                    if(allIntersects(nn,j,k) < initSphCoords(n,j,k) && ...
                            allIntersects(nn+1,j,k) > initSphCoords(n,j,k) )
                        bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                        %check to see if boundary
                        if(n ~=1)
                            if(bFlagSphere(n-1,j,k))
                                bFlagSphere(n-1,j,k) = bFlagSphere(n-1,j,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        if(n ~= nNodes+1)
                            if(bFlagSphere(n+1,j,k))
                                bFlagSphere(n+1,j,k) = bFlagSphere(n+1,j,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(j~=1)
                            if(bFlagSphere(n,j-1,k))
                                bFlagSphere(n,j-1,k) = bFlagSphere(n,j-1,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(j ~= nRayY)
                            if(bFlagSphere(n,j+1,k))
                                bFlagSphere(n,j+1,k) = bFlagSphere(n,j+1,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(k ~= 1)
                            if(bFlagSphere(n,j,k-1))
                                bFlagSphere(n,j,k-1) = bFlagSphere(n,j,k-1) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(k ~= nRayZ)
                            if(bFlagSphere(n,j,k+1))
                                bFlagSphere(n,j,k+1) = bFlagSphere(n,j,k+1) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                    end
                    
                end
            end
        end
        
    end
end

%Now need to reduce the number of spehres using a Greedy Algorithm
%but first find the boudary spheres and all internal spheres (including
%boudnary spheres)
%Also linearize the indices which makes your life easier later on
nbSph = sum(sum(sum(bFlagSphere >0 & bFlagSphere <7)));
nSph = sum(sum(sum(bFlagSphere >0)));
bndSph = zeros(3,nbSph);
intSph = zeros(3,nSph);
bndFlag= zeros(1,nSph);
%i2bIndexMap = zeros(1,nSph); %Map internal index to boundary index
tmp1 = 0;
tmp2 = 0;
for k = 1:nRayZ
    for j = 1:nRayY
        for n = 1:nNodes+1
            if(bFlagSphere(n,j,k)>0) 
                tmp1 = tmp1 + 1;
                intSph(:,tmp1) = [initSphCoords(n,j,k);initCoords(2,j,k); ...
                                  initCoords(3,j,k)];
                if(bFlagSphere(n,j,k)<7)
                    tmp2 = tmp2 + 1;
                    bndSph(:,tmp2) = [initSphCoords(n,j,k);initCoords(2,j,k); ...
                                      initCoords(3,j,k)];
                    bndFlag(tmp1) = 1;
                    %i2bIndexMap(tmp1)=tmp2;
                end
            end
        end
    end
end

if(tmp1 ~= nSph || tmp2 ~= nbSph)
    error('Something went wrong check the code!');
end

niSph = nSph - nbSph;
caSph = zeros(4,niSph);
tmp1  = 0;
for k = 1:nSph
    %Find coordinates of the candidate and its radius
    if(~bndFlag(k))
        tmp1 = tmp1 + 1;
        caCoords = intSph(:,k);
        
        dtmp = sqrt(sum(bsxfun(@minus,bndSph,caCoords).^2));
        %caRadius = min( setdiff(dtmp,min(dtmp)) )+ smoothFact*rSph;
        caRadius = min( dtmp ) + smoothFact*rSph;
        if(caRadius < rSph)
            warning('Candidate Radius is less than cell radius.');
        end
        %caRadius = dtmp(2)+eps(rSph);
        caSph(:,tmp1) = [caCoords; caRadius];
    end
    
end

if(tmp1 ~= niSph)
    error('Something went wrong check the code!');
end

nn = 0;
for k = 1:niSph
    nn = nn + sum((caSph(4,k)>sqrt(sum(bsxfun(@minus,intSph,caSph(1:3,k)).^2))));
end

sph2CandidateMap = spalloc(nSph, niSph, nn); 
for k = 1:niSph
    sph2CandidateMap(:,k) = (caSph(4,k)>sqrt(sum(bsxfun(@minus,intSph,caSph(1:3,k)).^2)))';
end

% for k = 1:niSph
%     for j = 1:nSph
%         origInds = sub2ind(size(sph2CandidateMap),j,k);
%     end
% end

%Array containing the original index (\in {1..niSph}) 
%of the final candidate sphere to be included in the cluster 
%finalClusterCellList array contains the list of indices of 
%original cells (\in {1..nSph}) exclusively contributing to any candidate 
%if there is an overlap the cells are added to the maximum candidate
%nCellPerSphere: counts the actual number of cells contributing to the candidate 
%i.e. (non-zero elements of the first dimension of finalClusterCellList)
%Assume initiall only 5.5% of the niSph will be include to the final list
finalClusterList = zeros(1,niSph); 
finalClusterCellList = zeros(nSph, fix(niSph/200));
nCellPerSphere = zeros(fix(niSph/200),1);


nclSph = 0; %number of final spheres from the niSph set to be included in the cluster

iinds = true(nSph,1);
jinds = true(niSph,1);

iindsOrg = linspace(1,nSph,nSph);
jindsOrg = linspace(1,niSph,niSph);


sz = sum(iinds);

while(sz > 0)
    
    %reduce all arrays with same sub indexing
    workMapArr = sph2CandidateMap(iinds,jinds);
    workIInds = iindsOrg(iinds);
    workJInds = jindsOrg(jinds);
    %Some tests to see if everything's Ok
    if(~sum(workMapArr(1,:))) 
        error('Spheres %d is assigned to non of the candidates.',k)
    end
    
    if(sum(workMapArr(1,:))==1)
        nclSph = nclSph + 1;
        j = find(workMapArr(1,:)==1);

        finalClusterList(nclSph)=workJInds(j);
         
         %remove the index of the current candidate for the next iteration
         jinds(workJInds(j)) = false;
         
         %remove the index of all cells within this candidate for the next
         %iteration
         ii = workMapArr(:,j) == 1;
         iinds(workIInds(ii)) = false;
        
        if(nclSph<=length(nCellPerSphere))
            nCellPerSphere(nclSph) = sum(ii);
            finalClusterCellList(1:nCellPerSphere(nclSph), nclSph) = workIInds(ii);
        else
            tmpArr1 = nCellPerSphere;
            tmpArr2 = finalClusterCellList;
            
            %Add another 50 elements
            finalClusterCellList = zeros(nSph, fix(niSph/200)+50);
            nCellPerSphere = zeros(fix(niSph/200)+50,1);
            
            nCellPerSphere(nclSph-1) = tmpArr1;
            finalClusterCellList(:, 1:nclSph-1) = tmpArr2;
            
            nCellPerSphere(nclSph) = sum(ii);
            finalClusterCellList(1:nCellPerSphere(nclSph), nclSph) = workIInds(ii);
        end
        
    else
        %find the max column
        [~,maxind] = max(sum(workMapArr,1));
        
        nclSph = nclSph + 1;
        %Add the original index to the current set of sphere cluster     
        finalClusterList(nclSph)=workJInds(maxind);
        
        %Remove the indices for the next iteration.
         jinds(workJInds(maxind)) = false;
         ii = workMapArr(:,maxind) == 1;
         iinds(workIInds(ii)) = false;    
        
        if(nclSph<=length(nCellPerSphere))
            nCellPerSphere(nclSph) = sum(ii);
            finalClusterCellList(1:nCellPerSphere(nclSph), nclSph) = workIInds(ii);
        else
            tmpArr1 = nCellPerSphere;
            tmpArr2 = finalClusterCellList;
            
            %Add another 50 elements
            sz = size(finalClusterCellList);
            finalClusterCellList = zeros(sz(1), sz(2)+50);
            nCellPerSphere = zeros(sz(2),1);
            
            nCellPerSphere(1:nclSph-1) = tmpArr1(1:nclSph-1);
            finalClusterCellList(:, 1:nclSph-1) = tmpArr2(:, 1:nclSph-1);
            
            nCellPerSphere(nclSph) = sum(ii);
            finalClusterCellList(1:nCellPerSphere(nclSph), nclSph) = workIInds(ii);
        end
            
    end
        
    sz = sum(iinds);  
        
end
 
assembly = caSph(:,finalClusterList(1:nclSph));
assembly(4,:)=2*assembly(4,:);
%now find the moment of inertia, contribution of individual spheres to the
%total mass and the volume of the assembly. 

%Assume cubic elements
nCellPerSphere = nCellPerSphere(1:nclSph);
cellMass = 8*rSph^3*rho; 
masses = cellMass*nCellPerSphere;
totMass = cellMass*sum(nCellPerSphere);

totVolSph = totMass/rho;

display(['STL Volume: ' num2str(stlVol) ' , Sphere Assembly Volume: ' num2str(totVolSph)]);
display(['Sum of candidate masses: ' num2str(sum(masses)) ' , total mass from all cells: ' num2str(totMass)]);

CoM = [0 0 0]';
for k = 1:nclSph
    inds = finalClusterCellList(1:nCellPerSphere(k), k);
    CoM = CoM + sum(intSph(1:3,inds),2);
end

%Multiply by cellMass and devide by totMass
CoM = CoM/sum(nCellPerSphere);

%Shift the CoM to origin
assembly(1:3,:) = bsxfun(@minus,assembly(1:3,:),CoM);
cellCentres = bsxfun(@minus,intSph,CoM);

CoMShifted = [0 0 0]';
for k = 1:nclSph
    inds = finalClusterCellList(1:nCellPerSphere(k), k);
    CoMShifted = CoMShifted + sum(cellCentres(1:3,inds),2);
end
CoMShifted = CoMShifted/sum(nCellPerSphere);

mI = zeros(3);
for k = 1:nclSph
    %Find the indices of the cells exclusively contributing to this cluster
    %sphere
    inds = finalClusterCellList(1:nCellPerSphere(k), k);
    coordMat = cellCentres(:,inds)*cellCentres(:,inds)';
    mI = mI + (trace(coordMat)*eye(3) - coordMat);
end

mI = cellMass*mI;

display('MI of the cluster: ');
display(mI);

if(findSTLmI)
    display('MI of the STL: ');
    display(STLmI);
end

display(['CoM of the cluster: ' num2str(CoM')])
if(findSTLmI)
    display(['CoM of the STL: ' num2str(STLCoM')]);
end
     
display(['CoM of the cluster after shifting: ' num2str(CoMShifted')]);

[filepath,fName,ext] = fileparts(fileName);
zBaseName = [fName '_SF' num2str(smoothFact)];
if(writeStlScaled)
    stlWrite([filepath fName 'Scaled' ext],faces, ...
                                           bsxfun(@minus,vertices,CoM'), ...
                                           'mode','ascii','title',fName{1})
end

if(writeTec)
    n = nRayZ*nRayY*(nNodes+1);
    zName = fullfile(filepath,[zBaseName '_BBoxSpheres' '.plt']);
    tecFid2=fopen(zName,'w');
    fprintf(tecFid2,'%s\n','Variables=');
    fprintf(tecFid2,'%s \n', ' "x" "y" "z" "d" ' );
    fprintf(tecFid2,'%s %s %s %d %s %d %s \n','ZONE T="',zBaseName, ...
            '" , STRANDID=',1,' ,I=',n, ...
            ' ZONETYPE=Ordered ,DATAPACKING=POINT');
    
    for k = 1:nRayZ
        for j = 1:nRayY
            for n = 1:nNodes+1
                fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f' '\n'],...
                    [initSphCoords(n,j,k)-CoM(1) ...
                    initCoords(2,j,k)-CoM(2) ...
                    initCoords(3,j,k)-CoM(3) rSph]);
                
                
            end
        end
    end
    
    fclose(tecFid2);
    zName = fullfile(filepath,[zBaseName '_STLDicrete' '.plt']);
    tecFid2=fopen(zName,'w');
    fprintf(tecFid2,'%s\n','Variables=');
    fprintf(tecFid2,'%s \n', ' "x" "y" "z" "d" "color"' );
    fprintf(tecFid2,'%s %s %s %d %s %d %s %d %s\n','ZONE T="',...
            zBaseName,'" , STRANDID=',1,' ,NODES=',8*nSph, ...
            ', ELEMENTS=' , nSph, ' ZONETYPE=FEBRICK ,DATAPACKING=POINT');
    
    for k = 1:nSph
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[-rSph;rSph;rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[-rSph;-rSph;rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[-rSph;-rSph;-rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[-rSph;rSph;-rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[rSph;rSph;rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[rSph;-rSph;rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[rSph;-rSph;-rSph];rSph;bndFlag(k)]');
        fprintf(tecFid2,['%16.10f %16.10f %16.10f %16.10f %d' '\n'],...
            [cellCentres(:,k)+[rSph;rSph;-rSph];rSph;bndFlag(k)]');
    end
    
    for k = 1:nSph
        fprintf(tecFid2,['%d %d %d %d %d %d %d %d' '\n'],...
            8*(k-1)+[1:8]');
    end
    
    fclose(tecFid2);
    
    zName = fullfile(filepath,[zBaseName '_STLDicrete' '.vtk']);
    tecFid2=fopen(zName,'w');
    fprintf(tecFid2,'# vtk DataFile Version 2.0\n');
    fprintf(tecFid2,[zBaseName '\n']);
    fprintf(tecFid2,'ASCII\n');
    fprintf(tecFid2,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(tecFid2,'\nPOINTS %d float\n', 8*nSph);
    
    for k = 1:nSph
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[-rSph;-rSph;-rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[rSph;-rSph;-rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[-rSph;rSph;-rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[rSph;rSph;-rSph]);
        
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[-rSph;-rSph;rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[rSph;-rSph;rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[-rSph;rSph;rSph]);
        fprintf(tecFid2,'%16.10e %16.10e %16.10e\n',...
                cellCentres(:,k)+[rSph;rSph;rSph]);
    end

    fprintf(tecFid2,'\nCELLS %d %d\n', nSph, 9*nSph);
    for k = 1:nSph
        fprintf(tecFid2,'%d %d %d %d %d %d %d %d %d\n',...
            [8;8*(k-1)+(1:8)'-1]);
    end    

    fprintf(tecFid2,'\nCELL_TYPES %d\n', nSph);
    fprintf(tecFid2,'%d\n',11*ones(1,nSph));
    
    fprintf(tecFid2,'\nCELL_DATA %d\n', nSph);
    fprintf(tecFid2,'\nSCALARS color int 1\n');
    fprintf(tecFid2,'LOOKUP_TABLE default\n');
    fprintf(tecFid2,'%d\n',bndFlag);
    
    fprintf(tecFid2,'\nSCALARS diam float 1\n');
    fprintf(tecFid2,'LOOKUP_TABLE default\n');
    fprintf(tecFid2,'%16.10e\n',rSph*ones(1,nSph));
    
    fclose(tecFid2);
    
    zName = fullfile(filepath,[zBaseName '_Cluster' '.plt']);
    tecFid2=fopen(zName,'w');
    fprintf(tecFid2,'%s\n','Variables=');
    fprintf(tecFid2,'%s \n', ' "x" "y" "z" "d" ' );
    fprintf(tecFid2,'%s %s %s %d %s %d %s \n','ZONE T="',zBaseName, ...
            '" , STRANDID=',1, ' ,I=',nclSph,...
            ' ZONETYPE=Ordered ,DATAPACKING=POINT');
    
    fprintf(tecFid2,'%16.10f %16.10f %16.10f %16.10f\n',assembly);
        
    fclose(tecFid2);
    
end

if(writeLammpsTemp)
    zName = fullfile(filepath,[zBaseName '_Cluster' '.lammps']);
    tecFid2=fopen(zName,'w');
    fprintf(tecFid2,'#%s\n',zBaseName);
    fprintf(tecFid2,'%d %s\n',nclSph,'atoms');
    fprintf(tecFid2,'%16.10e %16.10e %16.10e %s\n',CoMShifted','com');
    fprintf(tecFid2,[repmat('%16.10e ',1,6) ' %s\n'],...
            mI([1 5 9 2 3 6]),'inertia');
    fprintf(tecFid2,'%16.10e %s\n\n',totMass,'mass');
    
    sz = size(assembly);
    fprintf(tecFid2,'Coords\n\n');
    fprintf(tecFid2,'%d %16.10e %16.10e %16.10e\n',...
            [linspace(1,sz(2),sz(2));assembly(1:3,:)]);
    
    fprintf(tecFid2,'\nTypes\n\n');
    fprintf(tecFid2,'%d %d\n',[linspace(1,sz(2),sz(2)); ...
                               lammpsType*ones(1,sz(2))]);
    
    fprintf(tecFid2,'\nDiameters\n\n');
    fprintf(tecFid2,'%d %16.10e\n',[linspace(1,sz(2),sz(2));assembly(4,:)]);
    
    fprintf(tecFid2,'\nMasses\n\n');
    fprintf(tecFid2,'%d %16.10e\n',[linspace(1,sz(2),sz(2));masses']);
    
end
end


function  [STLmI,STLCoM] = calculateSTLmI(vertices,faces,stlVol,rho,nSamples)

xlow=min(vertices(:,1));
ylow=min(vertices(:,2));
zlow=min(vertices(:,3));

xhigh=max(vertices(:,1));
yhigh=max(vertices(:,2));
zhigh=max(vertices(:,3));

Lx = xhigh-xlow;
Ly = yhigh-ylow;
Lz = zhigh-zlow;

myPool=parpool(16);
nWorkers = myPool.NumWorkers;
nSamples = ceil(nSamples/nWorkers)*nWorkers;
myNumSamples = nSamples/nWorkers;
spmd
    [mySamples,myNTested]=workerSampleGenerator(vertices,faces,myNumSamples);
end

samples = zeros(3,nSamples);
nTested = 0;
for i = 1:nWorkers
    samples(:,(i-1)*myNumSamples+1:i*myNumSamples)=mySamples{i};
    nTested = nTested + myNTested{i};
end

delete(myPool);

estimatedVol = nSamples/nTested*Lx*Ly*Lz;

disp(['STL Volume Theoretical: ' num2str(stlVol) ' ,STL Volume MC: ' num2str(estimatedVol)]);
%Assume uniform mass distribution
STLCoM = mean(samples,2);

%Shift Samples to origin
samples = bsxfun(@minus,samples,STLCoM);
delm = stlVol*rho/nSamples;
coordMat = samples*samples';
STLmI = delm*(trace(coordMat)*eye(3) - coordMat);

end

function [mySamples,nTested]=workerSampleGenerator(vertices,faces,myNumSamples)

xlow=min(vertices(:,1));
ylow=min(vertices(:,2));
zlow=min(vertices(:,3));

xhigh=max(vertices(:,1));
yhigh=max(vertices(:,2));
zhigh=max(vertices(:,3));

Lx = xhigh-xlow;
Ly = yhigh-ylow;
Lz = zhigh-zlow;

%For each ray determine the intersection
sz = size(faces);
nface= sz(1);

ncSamples = 0;
nTested = 0;
samples = zeros(3, myNumSamples);
%acceptInds = false(1,nSample);
while ncSamples < myNumSamples
    cSamples = bsxfun(@plus,[xlow;ylow;zlow],bsxfun(@times,rand(3,myNumSamples),[Lx;Ly;Lz]));
    for n = 1:myNumSamples
        nTested = nTested + 1;
        direction = [1 0 0]';
        origin    = cSamples(:,n);
        ncIntersect = 0;
        for facet=1:nface
            vert = vertices(faces(facet,:),:);
            [flag, ~, ~, t] = rayTriangleIntersection(origin, direction, vert(1,:)', vert(2,:)', vert(3,:)', eps(Lx));
            if(flag && t>0)
                ncIntersect = ncIntersect+1;
            end
        end
        if(mod(ncIntersect,2)) %is inside
            ncSamples = ncSamples+1;
            samples(:,ncSamples) = cSamples(:,n);
            if(ncSamples==myNumSamples)
                break
            end
        end
    end
end

mySamples = samples;
end
