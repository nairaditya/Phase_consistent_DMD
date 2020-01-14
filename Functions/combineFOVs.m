function out = combineFOVs(U,X,Y)

% combines fiedls of views into a big one, doing weighted overlapse

% Ustruct is all the FOVs to combine in the format
% U(1).u(:,:,t) = data matrix 1
% U(1).v(:,:,t) = data matrix 1 (optional)
% U(1).w(:,:,t) = data matrix 1 (optional)
% U(1).X = meshgrid of x values for a
% U(1).X = dito for y

% X and Y are meshgrid of desired final combined FOV


% it goes like this:
% Determine overlap indication matrix, onese where any FOVs are
% overlapping, zeros else
% Make edge disdance matrix for each FOV, size is of full FO
% multiple by overlap indication matrix to zero non-overlap areas
% Divide matrix my its max value. Make non-overllaping areas of FOV == 1

% combined FOV will be (FOV1*W1 + FOV2*W2 + FOV3*W3...)/(W1 + W2 + w3...)



fieldNames = fields(U);
fieldNames(strcmp(fieldNames,'X')) = [];
fieldNames(strcmp(fieldNames,'Y')) = [];


overlap = zeros(size(X));

for i = 1:length(U)
   F{i} = griddedInterpolant(U(i).X.',U(i).Y.',U(i).(fieldNames{1})(:,:,1).','linear','none'); % interpoaltion object for FOV i
   
   % This matrix is one where FOV i is, zero elsewhere
   FOVind{i} = F{i}(X,Y);
   FOVind{i}(~isnan(FOVind{i})) = 1;
   FOVind{i}(FOVind{i}~=1) = 0;
   
   % include in overlap matrix
   overlap = overlap + FOVind{i};
   
   % edge distance matrix for FOV i
   edgeDist{i} = edgeDistance(FOVind{i});
   
end

% Map of where data exits
dataMap = overlap;
dataMap(overlap~=0) = 1;
dataMap(overlap==0) = NaN;

% make overlap equal one where overlaps occur, zero elsewhere
overlap(overlap == 1) = 0;
overlap(overlap>0) = 1;

denom = zeros(size(X));
for i = 1:length(U)
    % calculate weights matrix
    weights{i} = overlap.*edgeDist{i};
    w = weights{i};
    w(isnan(w)) = 0;
    weights{i} = weights{i}./max(w(:));
    % make weights = 1 where there is no overlap and FOV has data
    weights{i}(weights{i} == 0 & FOVind{i} == 1) = 1;
    denom = denom + weights{i};
end


fieldNames = fields(U);
fieldNames(strcmp(fieldNames,'X')) = [];
fieldNames(strcmp(fieldNames,'Y')) = [];

for i = 1:length(fieldNames) % iterate through data fields of U
    out.(fieldNames{i}) = zeros(size(X,1),size(X,2),size(U(1).(fieldNames{i}),3));
    
    for j = 1:size(U(1).(fieldNames{i}),3) % iterate through time
        for k = 1:length(U) % iterate through FOVs
            F{k}.Values = U(k).(fieldNames{i})(:,:,j).';
            data = F{k}(X,Y);
            data(isnan(data)) = 0;
            out.(fieldNames{i})(:,:,j) = out.(fieldNames{i})(:,:,j) + data.*weights{k};
        end
        out.(fieldNames{i})(:,:,j) = out.(fieldNames{i})(:,:,j).*dataMap./denom;
    end
end
