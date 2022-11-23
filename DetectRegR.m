function [ UpDnR ] = DetectRegR( X1, E1, Ref, R, RscaleUp )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        if (mean(E1(:,R)) >= RscaleUp*mean(X1(:,R)))&& (mean(E1(:,R)) >= RscaleUp*mean(Ref(:,R)))% If R Up-Regulated
            UpDnR = 1;      % R up-regulated through Pairing
        elseif (mean(E1(:,R)) <= (1/RscaleUp)*mean(X1(:,R)))&& (mean(E1(:,R)) <= (1/RscaleUp)*mean(Ref(:,R)))% If R Down-Regulated
            UpDnR = 2;      % R down regulated through Pairing
        else
            UpDnR = 0; % Pairing failed
        end
end

