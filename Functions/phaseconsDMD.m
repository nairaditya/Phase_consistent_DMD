function err = phaseconsDMD(comp,phiFOV,F,t,phiV,overlap,coords,mv)

%% Reconstruct FOVs with given phase shifts

% phiV = [0, phiV]; % phases for each FOV (origional and shift)
%comp = {'U','V','W'};
for j = 1:length(phiFOV)
    for k = 1:length(comp)
        for i = 1:length(t)
            % fov(j).(comp{k})(:,:,i) = real(phiFOV(j).(comp{k}).bkAk.*cos(wk*t(i)+phiFOV(j).(comp{k}).phi0+phiV(j)));
            fov(j).(comp{k})(:,:,i) = mv(j).(comp{k});
            
            for m = 1:length(phiFOV(j).(comp{k}))
                fov(j).(comp{k})(:,:,i) = fov(j).(comp{k})(:,:,i) + ...
                    real(phiFOV(j).(comp{k})(1).bkAk.*cos(phiFOV(j).(comp{k})(1).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(1).phi0));
            end
        end
    end
end

%% Compute errors
%comp = {'U','V','W'};
err = 0;
for i = 1:length(phiV) % this FOV gets interpolated onto FOV j overlap region
    for j = i+1:length(phiV)
        if any(overlap(i,j).m(:)) % only run if there actually is overlap
            % get the x,y values from the FOV that overlaps this one
            L = logical(overlap(j,i).m); % logical array size of FOV j indicating overlap with i
            x = coords(j).x(any(L'),any(L));
            y = coords(j).y(any(L'),any(L));
            m1 = size(x,1); n1 = size(x,2);
            for k = 1:length(comp)
                utemp1 = zeros(m1,n1,length(t));
                for h = 1:length(t)
                    F{i}.Values = fov(i).(comp{k})(:,:,h)'; % values of FOV i interpolated onto overlap in FOV j coordinates
                    utemp1(:,:,h) = F{i}(x',y')';
                end
                utemp2 = fov(j).(comp{k})(any(L'),any(L),:);
                err = err + sqrt(sum((utemp1(:)-utemp2(:)).^2));
            end
            
        end
    end
end

% if numel(err)>1
% send
