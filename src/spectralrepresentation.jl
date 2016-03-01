function acf(Z1)
    fZ1 = fft(Z1)
    pZ1 = abs2(fZ1)
    fftshift(real(ifft(pZ1)))
end
