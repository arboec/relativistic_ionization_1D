function term = calculateCoefficients(R,A0,c,taylorCoefficient,taylorCoefficient2) %p = pKinetic -- A, p0 = pCanonical -- A=0

a = -taylorCoefficient - A0*R;
b = -taylorCoefficient2;
a = a / (R*c);

squareRoot = sqrt(1 + a.*a);
squareRootOver = 1 ./ squareRoot;
denominator = 1.0 ./ (c*c * (1+a.*a).*(1+a.*a).*(1+a.*a).*(1+a.*a));% july -- removed 1/R*R from here

% 1
zeroTerm = squareRoot;
zeroTermCorrection = (1+a.*a).*(1+a.*a).*squareRoot.*denominator;
zeroTerm = zeroTerm - 0.5*1j*b.*zeroTermCorrection;

% 2
firstTerm = a .* squareRootOver;
firstTermCorrection = -3.0 * a .* (1 + a.*a).*squareRoot .* denominator;
firstTerm = firstTerm - 0.5*1j*b.*firstTermCorrection;

%3
secondTerm = c*c*zeroTermCorrection;
secondTermCorrection = -3.0 * (1-4*a.*a) .* squareRoot .* denominator;
secondTerm = secondTerm - 0.5*1j*b.*secondTermCorrection;

term = [zeroTerm; firstTerm; secondTerm];

denominator = 1.0 ./ ((squareRoot+1).*(squareRoot+1).*(squareRoot+1) .* squareRoot .* squareRoot .* squareRoot .* squareRoot .* squareRoot);
% 1
zeroTermDarwin = c*c ./ (squareRoot.*(squareRoot + 1));
zeroTermCorrectionDarwin = 6 * (a.*a.*(a.*a + (squareRoot + 1)) - 0.5*(squareRoot + 1)) .* denominator; %coorect, but slow
% 2
firstTermDarwin = - c*c * a.*(1+2*squareRoot) .* denominator .* squareRoot .* squareRoot .*(squareRoot+1);
%3
secondTermDarwin = c*c * zeroTermCorrectionDarwin;

termDarwin = [zeroTermDarwin; firstTermDarwin; secondTermDarwin];
term = [term; termDarwin];

end