#! /usr/bin/env python
# ^ this line tells the OS to use whatever the system python is when running the script, i.e. ./demonstrate_psd.py

import numpy as np


def compute_1d_psd(input_array, sample_spacing=1.0):

    # Ok, so the first thing that we need to do is define the Nyquist frequency.
    # This is the sample rate divided by two, or the highest frequency that can be resolved
    # in a discrete Fourier transform. The frequency content of the FFT is defined from - to + the nyquist frequency
    # Sometimes there are factors of 2pi in here, but I believe that numpy / scipy does not do that.

    # These are the defined freq for a discrete Fourier transform given the array size and sample spacing.
    # Here I've called fftshift to center the zero freq in the middle, which is more intuitive to me
    sample_freq = np.fft.fftshift(np.fft.fftfreq(array.shape[0], d=sample_spacing))

    # Now compute a discrete Fourier Transform
    array_fft = np.fft.fftshift(np.fft.fft(input_array))

    # The power spectrum is the square of the Fourier transform
    power_spectrum = np.real(array_fft * np.conj(array_fft))

    # or equivalently power_spectrum = np.abs(array_fft)**2

    # Now return the 1d power spectrum and the frequency sample points
    return sample_freq, power_spectrum


def compute_2d_psd(input_array, sample_spacing=1.0):

    # This works like the 1d case, with some small tweaks.
    # Now the frequency (this is easiest to think of as a k vector) has x and y components,
    # so create
    sample_freq_x = np.fft.fftshift(np.fft.fftfreq(array.shape[0], d=sample_spacing))
    sample_freq_y = np.fft.fftshift(np.fft.fftfreq(array.shape[1], d=sample_spacing))

    # Instead of a 1d array of freq, we now have a 2d array of freq denoted by k_x and k_y
    freq_xx, freq_yy = np.meshgrid(sample_freq_x, sample_freq_y)

    # Now we need a 2 dimensional FFT
    array_fft = np.fft.fftshift(np.fft.fft2(input_array))

    # Numpy array rules allow us to do this identically to the 1d case
    power_spectrum = np.real(array_fft * np.conj(array_fft))

    # Return 1d freq and squared Fourier transform points
    return np.sqrt(freq_xx ** 2 + freq_yy ** 2), power_spectrum


if __name__ == '__main__':

    # Compute a one dimensional PSD

    # Some general questions you can go for:
    # 1) Change the format of the input (for example, add some low frequency sine wave) and see what happens to the PSD
    # 2) Understand how the speed of the FFT depends on the size of the input array, particularly small number * (power of 2) versus a large prime
    # 3) Bin the output in frequency paying attention to numpy vectorization and array slicing

    array = np.random.normal(size=(100))
    freq, psd = compute_1d_psd(array, sample_spacing=0.1)
    print(freq, psd)

    # Compute a 2d psd

    # Play with the format of the input and try to get some intuition for how things map from position space -> harmonic space
    # You can even do crazy things like taking power spectra of image files.
    # Write something to bin azimuthally, i.e. into bins like 1 < kx**2 + ky**2 < 2 paying attention to vectorization and array slicing
    array = np.random.normal(size=(100, 100))
    freq, psd = compute_2d_psd(array, sample_spacing=0.1)
    print(freq, psd)
