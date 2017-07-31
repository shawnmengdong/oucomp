function mixture =  addMixture(components)
% This function creates a mixture structure

mixture = struct();
mixture.bip = zeroBIP(components); %modify later for non-zero BIP
mixture.components = components;

