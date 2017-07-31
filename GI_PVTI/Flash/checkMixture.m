function [mixture,zero_index] = checkMixture(mixture)

zero_index = find(~mixture.mole_fraction); %find 0 elements in mcomposition

if isempty(zero_index)==0
    mixture.components(zero_index)=[]; %delete the zero components
    mixture.bip.EOScons(:,zero_index)=[];
    mixture.bip.EOScons(zero_index,:)=[]; %delete the bip for the component
    mixture.bip.EOStdep(:,zero_index)=[];
    mixture.bip.EOStdep(zero_index,:)=[];
    mixture.mole_fraction(zero_index) = []; %delete the mole fraction
end

end