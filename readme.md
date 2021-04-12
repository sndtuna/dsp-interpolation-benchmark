## Overview
This is a small experiment to check the performance of interpolation methods for signal processing applications.

I was curious about methods that go beyond linear interpolation, but aren't restricted to a fixed ratio like FIR-filter up/downsampling. I didn't want to just copy paste magic constants from somewhere so I also wrote some small comments if you are curious how the polynomial interpolators are derived.

You can just look at the results from running the bench mark on the dev machine below or you can run the program yourself.

## Instructions
You need to have Rust installed to run following commands at the root directory of the repro:
* `cargo run` to run the benchmark and output the result in the same terminal.
* `cargo test` to (re-)generate the impulse response files.

## Frequency Spectrum
There is no evaluation of frequency spectrums, but you can use the impulse response files that are included in the repro and open them in a program with a spectrum analyzer like Audacity. To get the most accurate spectrum display with Audacity, set the size equal to the length of the files and use a rectangular window (The files have a length that is a power of two for this purpose).
The impulse response files have a high sampling rate, because they represent a use case of interpolating a 48000Hz sample rate source.

## Results
Output on the dev machine (Ryzen 2700X):
```
Using test length of 262144 samples...

warmup1:                      10.0 ns/sample.
warmup2:                      10.2 ns/sample.
warmup3:                      10.8 ns/sample.
linear:                        5.1 ns/sample.
cubic:                         9.7 ns/sample.
quintic:                      16.3 ns/sample.
quintic pure lagrange:        13.8 ns/sample.
truncated sinc:               18.4 ns/sample.
hann windowed sinc:           44.7 ns/sample.
truncated sinc(sin approx.):  14.8 ns/sample.
```

