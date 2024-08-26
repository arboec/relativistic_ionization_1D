function w = darwinTerm_formula(xi,RR,c,mu,firstStage,secondStage,thirdStage,p0,term,t,v,L) %p = pKinetic -- A, p0 = pCanonical -- A=0

p = p0;
p = p / (RR*c);

zeroTerm = term(4,:);
firstTerm = term(5,:);
secondTerm = term(6,:);

Vder2 = calcVder2(t,v,xi);
secondSlag = Vder2 .* (zeroTerm .*firstStage + mder1(xi,L).*secondStage.*firstTerm + 0.5*thirdStage .* (-1j*mder2(xi,L).*firstTerm /(RR*c) + mder1(xi,L).*mder1(xi,L).*secondTerm));
%===

firstStageP = fft(mu.*Vder2);
%firstStage = ifft(firstStageP);
firstStage = mu.*Vder2;

secondStageP = p .* firstStageP;
secondStage = ifft(secondStageP);

thirdStageP = p.*p .* firstStageP;
thirdStage = ifft(thirdStageP);

firstSlag = (zeroTerm .*firstStage + mder1(xi,L).*secondStage.*firstTerm + 0.5*thirdStage .* (-1j*mder2(xi,L).*firstTerm /(RR*c) + mder1(xi,L).*mder1(xi,L).*secondTerm));

w = 1/8*(firstSlag + secondSlag) / (c*c*c*c);

end
