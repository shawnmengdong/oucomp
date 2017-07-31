function BIP =  zeroBIP(components)

BIP = struct();
n = length(components);
BIP.EOScons = zeros(n);
BIP.EOStdep = zeros(n);