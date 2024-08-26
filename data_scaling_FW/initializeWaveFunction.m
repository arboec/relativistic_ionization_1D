function mu = initializeWaveFunction(FILENAME_RE,FILENAME_IM,N)

    fileName = FILENAME_RE;
    fid = fopen(fileName, 'r');
    muRe = fscanf(fid, '%f ', N);
    fclose(fid);
    fileName = FILENAME_IM;
    fid = fopen(fileName, 'r');
    muIm = fscanf(fid, '%f ', N);
    fclose(fid);  
    mu = muRe + 1j*muIm;

end