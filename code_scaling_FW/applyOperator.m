function w = applyOperator(t,v,xi,c,mu,p0Xi,term,L) %p = pKinetic -- A, p0 = pCanonical -- A=0

RR = R(t,v);
p = p0Xi / (RR*c);

%part with mu
p = fftshift(p);
muP = fft(mu);

secondStageP = p .* muP;
secondStage = ifft(secondStageP);

thirdStageP = p.*p .* muP;
thirdStage = ifft(thirdStageP);

w = -1j*c*c * (mder1(xi,L).*secondStage.*term(2,:) + 0.5*term(3,:) .* (-1j*mder2(xi,L).*secondStage/(RR*c) + mder1(xi,L).*mder1(xi,L).*thirdStage));
w = w - 0*1j*darwinTerm_formula(xi,RR,c,mu,mu,secondStage,thirdStage,p0Xi,term,t,v,L); % darwinTerm
w = w + c*1j*Rd(t,v)*xi .* mder1(xi,L).*secondStage;
w = -w;
end