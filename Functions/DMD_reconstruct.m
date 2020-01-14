function fov = DMD_reconstruct(comp,phiFOV,phiV,mv,t)

%comp = {'U','V','W'};
for j = 1:length(phiFOV)
    for k = 1:length(comp)
        for i = 1:length(t)
            if length(comp) == 3
            if ~isempty(mv)
                % fov(j).(comp{k})(:,:,i) = real(phiFOV(j).(comp{k}).bkAk.*cos(wk*t(i)+phiFOV(j).(comp{k}).phi0+phiV(j)));
                fov(j).(comp{k})(:,:,i) = mv(j).(comp{k}) + real(phiFOV(j).(comp{k})(1).bkAk.*cos(phiFOV(j).(comp{k})(1).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(1).phi0)+...
                    phiFOV(j).U(2).bkAk.*cos(phiFOV(j).(comp{k})(2).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(2).phi0));
            else
                fov(j).(comp{k})(:,:,i) = real(phiFOV(j).(comp{k})(1).bkAk.*cos(phiFOV(j).(comp{k})(1).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(1).phi0)+...
                    phiFOV(j).U(2).bkAk.*cos(phiFOV(j).(comp{k})(2).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(2).phi0));
            end
            elseif length(comp) == 1
                if ~isempty(mv)
                % fov(j).(comp{k})(:,:,i) = real(phiFOV(j).(comp{k}).bkAk.*cos(wk*t(i)+phiFOV(j).(comp{k}).phi0+phiV(j)));
                fov(j).(comp{k})(:,:,i) = mv(j).(comp{k}) + real(phiFOV(j).(comp{k})(1).bkAk.*cos(phiFOV(j).(comp{k})(1).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(1).phi0)+...
                    phiFOV(j).VORT(2).bkAk.*cos(phiFOV(j).(comp{k})(2).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(2).phi0));
            else
                fov(j).(comp{k})(:,:,i) = real(phiFOV(j).(comp{k})(1).bkAk.*cos(phiFOV(j).(comp{k})(1).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(1).phi0)+...
                    phiFOV(j).VORT(2).bkAk.*cos(phiFOV(j).(comp{k})(2).wk*(t(i)+phiV(j)) + phiFOV(j).(comp{k})(2).phi0));
                end
            end
        end
    end
end


% function st = DMD_reconstruct(Phi,omega,b,phiOff,m,n,t)
%
% eind = 0;
%
%      for j = 1:3 % through components
%         sind = eind + 1;
%         eind = eind + m(1)*n(1);
%         ind1(j) = sind;
%         ind2(j) = eind;
%     end
%
% cmp = {'U','V','W'};
%
%
% for i = 1:length(m) % through FOVs
%
%     for j = 1:3 % through components
%         st(i).(cmp{j}) = zeros(m(i),n(i),length(t));
%     end
%
%
%     for k = 15:16%:length(b) % iterate mode
%
%         a = real(Phi(:,k));
%         c = imag(Phi(:,k));
%         A = abs(Phi(:,k)); % Amplitude
%
%         bkAk = b(k).*A;
%         wk = imag(omega(k));
%         phi0 = 2*atan2(c,a+A);
%
%
%         for p = 1:length(t) % iterate time, calculating mode-specific X, with FOV-specific phase offset
%             % XTemp(:,p) = real(b(k).*A.*cos(wk*t(p)+phi0+phiOff(k,i)));
%             XTemp(:,p) = real(b(k).*A.*cos(wk*(t(p)+phiOff(i))+phi0));
%         end
%
%         % add to the components
%
%         for j = 1:3 % through components
%
%             st(i).(cmp{j}) = st(i).(cmp{j}) + reshape(XTemp(ind1(j):ind2(j),:),m(i),n(i),length(t));
%         end
%
%
%
%     end
%
%      for j = 1:3 % through components
%         sind = eind + 1;
%         eind = eind + m(i)*n(i);
%         ind1(j) = sind;
%         ind2(j) = eind;
%     end
%
%
% end








end